export Compute_Corrections_Init, Compute_Corrections!, Four_In_One!, Spectral_Dynamics!, Get_Topography!, Spectral_Initialize_Fields!, Spectral_Dynamics_Physics!, Atmosphere_Update!

function Compute_Corrections_Init(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
                                  grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3}, grid_ps_p::Array{Float64, 3}, grid_t_p::Array{Float64, 3}, 
                                  grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_δt::Array{Float64, 3},  
                                  Δt::Int64, grid_energy_temp::Array{Float64, 3}, grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_δtracers::Array{Float64,3})
    """
    Computes baseline global mass and energy integrals for conservation enforcement.

    This function serves as the initialization step for the model's a posteriori conservation 
    schemes. It calculates the reference global states (mass and total energy) that the model 
    attempts to preserve or restore during the correction phase (Compute_Corrections!).

    It performs the following operations:
    1.  Mass Reference: Computes the area-weighted global mean surface pressure (`mean_ps_p`) 
        if mass correction is enabled.
    2.  Energy Reference: Computes the mass-weighted global integral of total energy 
        (kinetic + enthalpy) if energy correction is enabled. 
        Note: The energy calculation projects the state forward using the provided tendencies 
        (grid_δ variables) and time step (Δt) to establish the target energy level.
    3.  Validation: Checks for unsupported correction modes (e.g., water correction) and 
        raises errors if triggered.

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        The vertical coordinate infrastructure defining layers and integration weights.

    mesh (Spectral_Spherical_Mesh)
        The spherical mesh geometry containing area weights and grid dimensions.

    atmo_data (Atmo_Data)
        Structure containing physical constants (cp_air, grav) and configuration flags 
        (do_mass_correction, do_energy_correction).

    grid_u_p (Array{Float64, 3})
        Zonal wind component at the previous time step (or predictor state).

    grid_v_p (Array{Float64, 3})
        Meridional wind component at the previous time step (or predictor state).

    grid_ps_p (Array{Float64, 3})
        Surface pressure field at the previous time step.

    grid_t_p (Array{Float64, 3})
        Temperature field at the previous time step.

    grid_δu (Array{Float64, 3})
        Tendency of the zonal wind component.

    grid_δv (Array{Float64, 3})
        Tendency of the meridional wind component.

    grid_δt (Array{Float64, 3})
        Tendency of the temperature field.

    Δt (Int64)
        The integration time step in seconds.

    grid_energy_temp (Array{Float64, 3})
        A pre-allocated buffer array used to store the intermediate energy field calculation 
        to avoid memory allocation.

    Returns
    -------
    mean_ps_p (Float64)
        The global area-weighted mean surface pressure (Pascal).

    mean_energy_p (Float64)
        The global mass-weighted integral of total energy (Joules/m^2).

    Modified
    --------
    grid_energy_temp
        This array is overwritten with the calculated specific total energy field 
        (Kinetic + Enthalpy) during execution.

    """
    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    if (do_mass_correction) 
        mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
    end
    
    if (do_energy_correction) 
        # due to dissipation introduced by the forcing
        cp_air, grav       = atmo_data.cp_air, atmo_data.grav
        grid_energy_temp  .= 0.5*((grid_u_p + Δt*grid_δu).^2 + (grid_v_p + Δt*grid_δv).^2) + cp_air*(grid_t_p + Δt*grid_δt)
        mean_energy_p      = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_p)
    end
    
    if (do_water_correction)
        # error("water correction has not implemented")
        mean_moisture_p    =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_p .+ grid_δtracers * Δt, grid_ps_p)
    end
    
    return mean_ps_p, mean_energy_p, mean_moisture_p 
end 

function Compute_Corrections!(semi_implicit::Semi_Implicit_Solver, vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
                              mean_ps_p::Float64, mean_energy_p::Float64, mean_moisture_p::Float64,
                              grid_u_n::Array{Float64, 3}, grid_v_n::Array{Float64, 3},
                              grid_energy_temp::Array{Float64, 3}, grid_ps_p::Array{Float64, 3},grid_ps_c::Array{Float64, 3},
                              grid_ps_n::Array{Float64, 3}, spe_lnps_n::Array{ComplexF64, 3}, 
                              grid_t_n::Array{Float64, 3}, spe_t_n::Array{ComplexF64, 3},
                              grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_tracers_n::Array{Float64, 3}, 
                              grid_t::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, grid_z_full::Array{Float64, 3}, grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3},
                              grid_geopots::Array{Float64, 3}, grid_w_full::Array{Float64,3}, grid_t_p::Array{Float64, 3}, dyn_data::Dyn_Data, grid_δt::Array{Float64,3})#,  grid_tracers_diff::Array{Float64, 3})
    """
    Applies global conservation fixes to the updated atmospheric state variables.

    This function acts as the corrector step in the model's conservation scheme. It modifies 
    the prognostic variables for the next time step (denoted by suffix `_n`) to ensure that 
    the global mass and total energy match the reference values calculated previously in 
    `Compute_Corrections_Init`.

    It executes the following correction algorithms:
    1.  Mass Correction:
        - Calculates the global mean surface pressure of the updated field (`grid_ps_n`).
        - Derives a multiplicative scaling factor (`mass_correction_factor = target / current`).
        - Scales the entire grid-point surface pressure field.
        - Updates the (0,0) spectral coefficient of log-pressure (`spe_lnps_n`) to reflect 
        this shift.
    2.  Energy Correction:
        - Computes the global total energy of the updated state using new winds (`grid_u_n`, 
        `grid_v_n`) and temperatures.
        - Calculates the global energy deficit/surplus.
        - Derives a uniform temperature increment (`temperature_correction`) required to close 
        the energy budget.
        - Adds this increment to the grid-point temperature (`grid_t_n`) and the zonal mean 
        spectral temperature coefficients (`spe_t_n`).

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        The vertical coordinate infrastructure for integration weights.

    mesh (Spectral_Spherical_Mesh)
        The spherical mesh geometry containing area weights.

    atmo_data (Atmo_Data)
        Structure containing physical constants (cp_air, grav) and configuration flags.

    mean_ps_p (Float64)
        The target global mean surface pressure (conserved quantity from previous step).

    mean_energy_p (Float64)
        The target global total energy integral (conserved quantity from previous step).

    grid_u_n (Array{Float64, 3})
        Updated zonal wind component (time t + Δt).

    grid_v_n (Array{Float64, 3})
        Updated meridional wind component (time t + Δt).

    grid_energy_temp (Array{Float64, 3})
        A pre-allocated buffer array used for intermediate energy calculations.

    grid_ps_n (Array{Float64, 3})
        Updated surface pressure field (time t + Δt).

    spe_lnps_n (Array{ComplexF64, 3})
        Updated spectral coefficients of log-surface pressure (time t + Δt).

    grid_t_n (Array{Float64, 3})
        Updated temperature field (time t + Δt).

    spe_t_n (Array{ComplexF64, 3})
        Updated spectral coefficients of temperature (time t + Δt).

    Returns
    -------
    None
        The function operates in-place on the `_n` arrays.

    Modified
    --------
    grid_ps_n
        Scaled in-place to enforce mass conservation.

    spe_lnps_n
        The global mean component (0,0) is adjusted to match the pressure scaling.

    grid_t_n
        Shifted uniformly in-place to enforce energy conservation.

    spe_t_n
        The zonal mean components are shifted to match the temperature correction.

    grid_energy_temp
        Overwritten with the energy field of the uncorrected state.

    """
    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    
    if (do_mass_correction) 
        mean_ps_n              = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        mass_correction_factor = mean_ps_p/mean_ps_n
        grid_ps_n            .*= mass_correction_factor
        #P00 = 1 
        spe_lnps_n[1,1,1]     += log(mass_correction_factor)
    end
    
    if (do_energy_correction) 
        cp_air, grav           = atmo_data.cp_air, atmo_data.grav
        grid_energy_temp      .= 0.5*(grid_u_n.^2 + grid_v_n.^2) + cp_air*grid_t_n
        mean_energy_n          = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_n)
        
        temperature_correction = grav*(mean_energy_p - mean_energy_n)/(cp_air*mean_ps_p)
        #@info grav, mean_energy_p , mean_energy_n, cp_air, mean_ps_p
        grid_t_n             .+= temperature_correction
        spe_t_n[1,1,:]       .+= temperature_correction
    end

    # @info mean_ps_p, mean_energy_p, mass_correction_factor, temperature_correction
    # error(6868)
    
    if (do_water_correction) 
        ### By CJY 0517
        grid_tracers_n[grid_tracers_n .< 0.] .=  0.
        mean_moisture_n                       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        grid_tracers_n                      .*=  mean_moisture_p ./ mean_moisture_n 
        mean_moisture_n                       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)        
        ### 10/30 
        # @info "#### moisture correction:", (mean_moisture_n - mean_moisture_p)
        return mean_moisture_n
    end
  
end 



"""
compute vertical mass flux and velocity 
grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
grid_w_full[:,:,k]   = dp/dt vertical velocity 

update residuals
grid_δps[:,:,k]  += -∑_{r=1}^nd Dr = -∑_{r=1}^nd ∇(vrΔp_r)
grid_δt[:,:,k]  += κTw/p 
(grid_δu[:,:,k], grid_δv[:,:,k]) -= RT ∇p/p 

!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: grid_M_half = grid_ps*d(sigma)/dt

"""
function Four_In_One!(vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
                      grid_div::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
                      grid_ps::Array{Float64,3},  grid_Δp::Array{Float64,3}, grid_lnp_half::Array{Float64,3}, grid_lnp_full::Array{Float64,3}, grid_p_full::Array{Float64,3},
                      grid_dλ_ps::Array{Float64,3}, grid_dθ_ps::Array{Float64,3}, 
                      grid_t::Array{Float64,3}, 
                      grid_M_half::Array{Float64,3}, grid_w_full::Array{Float64,3}, 
                      grid_δu::Array{Float64,3}, grid_δv::Array{Float64,3}, grid_δps::Array{Float64,3}, grid_δt::Array{Float64,3}, grid_δtracers::Array{Float64,3})
    """
    Computes diagnostic vertical kinematic variables and updates prognostic tendencies with vertical coupling terms.

    This function aggregates four distinct but structurally related dynamical calculations into a 
    single vertical loop to maximize computational efficiency. It implements the generalized 
    vertical discretization scheme proposed by Simmons and Burridge (1981) for hybrid sigma-pressure 
    coordinates.

    It performs the following "Four" key operations:
    1.  Vertical Mass Flux (M): Computes the downward mass flux per unit area across layer interfaces 
        required for vertical advection.
    2.  Vertical Velocity (ω): Calculates the pressure vertical velocity (Dp/Dt) at full model levels.
    3.  Pressure Gradient Force: Computes the solenoidal term (-RT ∇lnp) and adds it to the 
        momentum tendencies (grid_δu, grid_δv).
    4.  Energy Conversion: Computes the adiabatic expansion/compression term (κTω/p) and adds it 
        to the temperature tendency (grid_δt).

    Additionally, it integrates the column divergence to update the surface pressure tendency (grid_δps).

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        Vertical grid definition containing coefficients (ak, bk) and difference options.

    atmo_data (Atmo_Data)
        Physical constants including Gas constant (rdgas) and specific heat (cp_air).

    grid_div (Array{Float64,3})
        Divergence of the horizontal wind field on the grid.

    grid_u (Array{Float64,3})
        Zonal wind component.

    grid_v (Array{Float64,3})
        Meridional wind component.

    grid_ps (Array{Float64,3})
        Surface pressure.

    grid_Δp (Array{Float64,3})
        Pressure thickness of each model layer.

    grid_lnp_half (Array{Float64,3})
        Logarithm of pressure at layer interfaces.

    grid_lnp_full (Array{Float64,3})
        Logarithm of pressure at full layer midpoints.

    grid_p_full (Array{Float64,3})
        Pressure at full layer midpoints.

    grid_dλ_ps (Array{Float64,3})
        Zonal gradient of surface pressure (∂ps/∂λ).

    grid_dθ_ps (Array{Float64,3})
        Meridional gradient of surface pressure (∂ps/∂θ).

    grid_t (Array{Float64,3})
        Temperature field.

    grid_M_half (Array{Float64,3})
        Output buffer for vertical mass flux at interfaces.

    grid_w_full (Array{Float64,3})
        Output buffer for vertical pressure velocity at full levels.

    grid_δu (Array{Float64,3})
        Accumulator for zonal wind tendency.

    grid_δv (Array{Float64,3})
        Accumulator for meridional wind tendency.

    grid_δps (Array{Float64,3})
        Accumulator for surface pressure tendency.

    grid_δt (Array{Float64,3})
        Accumulator for temperature tendency.

    Returns
    -------
    None
        Operates in-place on the output and tendency arrays.

    Modified
    --------
    grid_M_half
        Filled with the calculated vertical mass flux.

    grid_w_full
        Filled with the calculated vertical velocity (ω).

    grid_δu, grid_δv
        Decremented by the pressure gradient force (-RT ∇lnp).

    grid_δt
        Incremented by the energy conversion term (κTω/p).

    grid_δps
        Decremented by the integrated column divergence (-∑∇(vΔp)).

    """
    rdgas, cp_air          = atmo_data.rdgas, atmo_data.cp_air
    nd, bk                 = vert_coord.nd, vert_coord.bk
    Δak, Δbk               = vert_coord.Δak, vert_coord.Δbk
    vert_difference_option = vert_coord.vert_difference_option
    
    kappa                  = rdgas / cp_air
    
    # dmean_tot = ∇ ∑_{k=1}^{nd} vk Δp_k = ∑_{k=1}^{nd} Dk
    nλ, nθ, _              = size(grid_ps)
    dmean_tot              = zeros(Float64, nλ, nθ)
    Δlnp_p                 = zeros(Float64, nλ, nθ)
    Δlnp_m                 = zeros(Float64, nλ, nθ)
    Δlnp                   = zeros(Float64, nλ, nθ)
    x1                     = zeros(Float64, nλ, nθ)
    dlnp_dλ                = zeros(Float64, nλ, nθ)
    dlnp_dθ                = zeros(Float64, nλ, nθ)
    dmean                  = zeros(Float64, nλ, nθ)
    x5                     = zeros(Float64, nλ, nθ)
      
    if (vert_difference_option == "simmons_and_burridge") 
        for k = 1:nd
            Δp       = grid_Δp[:,:,k]
            
            Δlnp_p  .= grid_lnp_half[:,:,k + 1] - grid_lnp_full[:,:,k]
            Δlnp_m  .= grid_lnp_full[:,:,k]     - grid_lnp_half[:,:,k]
            Δlnp    .= grid_lnp_half[:,:,k + 1] - grid_lnp_half[:,:,k]
            
            # angular momentum conservation 
            #    ∇p_k/p =  [(lnp_k - lnp_{k-1/2})∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k)∇p_{k+1/2}]/Δpk
            #         =  [(lnp_k - lnp_{k-1/2})B_{k-1/2} + (lnp_{k+1/2} - lnp_k)B_{k+1/2}]/Δpk * ∇ps
            #         =  x1 * ∇ps
            x1      .= (bk[k] * Δlnp_m + bk[k + 1] * Δlnp_p ) ./ Δp
            
            dlnp_dλ .= x1 .* grid_dλ_ps[:,:,1]
            dlnp_dθ .= x1 .* grid_dθ_ps[:,:,1]
            
            
            
            # (grid_δu, grid_δv) -= RT ∇p/p 
            grid_δu[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dλ
            grid_δv[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dθ
            
            # dmean = ∇ (vk Δp_k) =  divk Δp_k + vk  Δbk[k] ∇ p_s
            dmean           .= grid_div[:,:,k] .* Δp + Δbk[k] * (grid_u[:,:,k] .* grid_dλ_ps[:,:,1] + grid_v[:,:,k] .* grid_dθ_ps[:,:,1])
            
        
            # energy conservation for temperature
            # w/p = dlnp/dt = ∂lnp/∂t + dσ ∂lnp/∂σ + v∇lnp
            # dσ ∂ξ_k/∂σ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k
            # weight the same way (TODO)
            # vertical advection operator (M is the downward speed)
            # dσ ∂lnp_k/∂σ = [M_{k+1/2}(lnp_k+1/2 - lnp_k) + M_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
            # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2}(lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
            # As we know
            # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
            
            # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
            #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
            
            x5                     .= -(dmean_tot .* Δlnp + dmean .* Δlnp_p) ./ Δp .+ grid_u[:,:,k] .* dlnp_dλ + grid_v[:,:,k] .* dlnp_dθ
            # grid_δt += κT w/p
            grid_δt[:,:,k]        .+=  kappa * grid_t[:,:,k] .* x5
            # grid_w_full = w
            grid_w_full[:,:,k]     .= x5 .* grid_p_full[:,:,k]
            # update dmean_tot to ∑_{r=1}^k ∇(vrΔp_r)
            dmean_tot             .+= dmean
            # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
            grid_M_half[:,:,k + 1] .= -dmean_tot
        end
      
    else
        error("vert_difference_option ", vert_difference_option, " is not a valid value for option")
      
    end
    # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
    grid_δps[:,:,1]        .-= dmean_tot
    
    for k = 1:nd-1
        # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
        grid_M_half[:,:,k+1] .+= dmean_tot * bk[k+1]
    end
    
    grid_M_half[:,:,1]      .= 0.0
    grid_M_half[:,:,nd + 1] .= 0.0

end 



"""
The governing equations are
∂div/∂t = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t= (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t = -(u,v)∇T - dσ∂T∂σ + κTw/p + J:= f^t           
Φ = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 

We have 
δdiv = f^d - I^d + I^d
δlnps = f^p - I^p + I^p
δT = f^t - I^t + I^t

"""
function Spectral_Dynamics!(mesh::Spectral_Spherical_Mesh,  vert_coord::Vert_Coordinate, 
                            atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
                            semi_implicit::Semi_Implicit_Solver)
    """
    Orchestrates the complete semi-implicit spectral time-stepping cycle for the atmospheric primitive equations.

    This function serves as the central dynamical core driver. It advances the atmospheric state 
    from time t to t + Δt by integrating the vorticity, divergence, temperature, and surface 
    pressure equations. It employs a hybrid approach, calculating non-linear advection and 
    physical forcings in grid space, while solving linear wave propagation and derivatives in 
    spectral space.

    The execution pipeline proceeds as follows:
    1.  Conservation Init: Captures the reference mass and energy states for later correction.
    2.  Grid-Point Dynamics:
        -   Computes 3D pressure fields and geopotential height.
        -   Calculates vertical velocities and thermodynamic terms via `Four_In_One!`.
        -   Computes vertical advection (finite volume/difference) and horizontal advection.
        -   Applies the Coriolis force and non-linear rotational terms to wind tendencies.
    3.  Spectral Transformation:
        -   Transforms grid-space tendencies (u, v, T, lnps) into spectral coefficients 
            (vorticity, divergence, T, lnps).
        -   Computes the Laplacian of total energy (Φ + K) required for the divergence equation.
    4.  Semi-Implicit Solver:
        -   Solves the coupled linear system for gravity waves and acoustic modes implicitly 
            to permit large time steps (`Implicit_Correction!`).
    5.  Time Integration:
        -   Applies horizontal diffusion (spectral damping).
        -   Advances the state variables using a Time-Filtered Leapfrog scheme (`Filtered_Leapfrog!`).
    6.  Recovery & Correction:
        -   Transforms the new spectral states back to grid space.
        -   Applies a posteriori mass and energy corrections (`Compute_Corrections!`).
        -   Updates internal time counters.

    Parameters
    ----------
    mesh (Spectral_Spherical_Mesh)
        Geometry and spectral transform infrastructure (Legendre polynomials, FFTs).

    vert_coord (Vert_Coordinate)
        Vertical grid definition.

    atmo_data (Atmo_Data)
        Physical constants and model configuration.

    dyn_data (Dyn_Data)
        Container for all prognostic variables (current `_c`, previous `_p`, next `_n`) 
        and intermediate work arrays.

    semi_implicit (Semi_Implicit_Solver)
        Structure containing the reference state profiles and matrices for the implicit 
        linear solver.

    Returns
    -------
    None
        The function advances the state of `dyn_data` in-place.

    Modified
    --------
    dyn_data.spe_X_n (X = vor, div, t, lnps)
        Updated spectral coefficients for the next time step.

    dyn_data.grid_X_n
        Updated grid-point fields for the next time step.

    dyn_data.grid_δX
        Tendency arrays are overwritten during the computation of non-linear terms.

    integrator.time
        Increments the simulation clock by Δt.

    """
    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n


    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps                                 = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps                                              = dyn_data.grid_lnps
    
    grid_div, grid_absvor, grid_vor                        = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half                               = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half       = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
    
    grid_energy_full, spe_energy                           = dyn_data.grid_energy_full, dyn_data.spe_energy
    
    # By CJY2
    spe_tracers_n     = dyn_data.spe_tracers_n
    spe_tracers_c     = dyn_data.spe_tracers_c
    spe_tracers_p     = dyn_data.spe_tracers_p 
    
    grid_tracers_n    = dyn_data.grid_tracers_n
    grid_tracers_c    = dyn_data.grid_tracers_c
    grid_tracers_p    = dyn_data.grid_tracers_p 
    
    grid_tracers_diff = dyn_data.grid_tracers_diff
    
    spe_δtracers      = dyn_data.spe_δtracers  
    grid_δtracers     = dyn_data.grid_δtracers 
    
    ### 11/07
    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    ###
    grid_w_full       = dyn_data.grid_w_full
    grid_δtracers     = dyn_data.grid_δtracers 
    ###
    # original 
    # todo !!!!!!!!
    #  grid_q = grid_t
    
    # pressure difference
    grid_Δp = dyn_data.grid_Δp
    # temporary variables
    grid_δQ = dyn_data.grid_d_full1
    
    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
    

    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)

    mean_ps_p, mean_energy_p, mean_moisture_p = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
    grid_u_p, grid_v_p, grid_ps_p, grid_t_p, 
    grid_δu, grid_δv, grid_δt,  
    Δt, grid_energy_full, grid_tracers_p, grid_tracers_c, grid_δtracers)
    
    # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    
    # compute ∇ps = ∇lnps * ps
    Compute_Gradients!(mesh, spe_lnps_c,  grid_dλ_ps, grid_dθ_ps)
    grid_dλ_ps .*= grid_ps
    grid_dθ_ps .*= grid_ps


    
    # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
    # except the contributions from geopotential or vertical advection
    Four_In_One!(vert_coord, atmo_data, grid_div, grid_u, grid_v, grid_ps, 
    grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
    grid_dλ_ps, grid_dθ_ps, 
    grid_t, 
    grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, grid_δtracers)
    
    Compute_Geopotential!(vert_coord, atmo_data, 
    grid_lnp_half, grid_lnp_full,  
    grid_t, 
    grid_geopots, grid_geopot_full, grid_geopot_half, grid_tracers_c)
    

    grid_δlnps .= grid_δps ./ grid_ps
    Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)
    

    
    # compute vertical advection, todo  finite volume method 
    Vert_Advection!(vert_coord, grid_u, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δu  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_v, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δv  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_t, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δt  .+= grid_δQ
    
    ### By CJY2 spectral tracers need to be done first 
    Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme,  grid_δQ)
    grid_δtracers .+= grid_δQ 
    Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u, grid_v, grid_δtracers) 
    Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
    Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
    Filtered_Leapfrog!(integrator, spe_δtracers, spe_tracers_p, spe_tracers_c, spe_tracers_n)
    Trans_Spherical_To_Grid!(mesh, spe_tracers_n, grid_tracers_n)
    ###################################################
  
    Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)
    Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)

    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    
    
    grid_δu .+=  grid_absvor .* grid_v
    grid_δv .-=  grid_absvor .* grid_u
    
    
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)
    

    grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u.^2 + grid_v.^2)
    Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy
    
    
    Implicit_Correction!(semi_implicit, vert_coord, atmo_data,
    spe_div_c, spe_div_p, spe_lnps_c, spe_lnps_p, spe_t_c, spe_t_p, 
    spe_δdiv, spe_δlnps, spe_δt)


    
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)

    
    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
    Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)
    

    
    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
    Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
    grid_ps_n .= exp.(grid_lnps)


    mean_moisture_n_loc = Compute_Corrections!(semi_implicit, vert_coord, mesh, atmo_data, mean_ps_p, mean_energy_p,mean_moisture_p, 
    grid_u_n, grid_v_n,
    grid_energy_full, grid_ps_p,grid_ps,
    grid_ps_n, spe_lnps_n, 
    grid_t_n, spe_t_n, 
    grid_tracers_p, grid_tracers_c, grid_tracers_n,
    grid_t, grid_p_full, grid_p_half, grid_z_full, grid_u_p, grid_v_p, grid_geopots, grid_w_full, grid_t_p, dyn_data, grid_δt)

   
    day_to_sec = 86400
    if (integrator.time%(day_to_sec/4) == 0)
        # dyn_data.grid_tracers_c[dyn_data.grid_tracers_c .< 0] .= 0
        @info "Day: ", (integrator.time/ day_to_sec), " Max |U|,|V|,|P|,|T|,|qv|: ", maximum(abs.(dyn_data.grid_u_c)), maximum(abs.(dyn_data.grid_v_c)), maximum(dyn_data.grid_p_full), maximum(dyn_data.grid_t_c), maximum(dyn_data.grid_tracers_c)
        @info "Day: ", (integrator.time/ day_to_sec), " Min |U|,|V|,|P|,|T|,|qv|: ", minimum(abs.(dyn_data.grid_u_c)), minimum(abs.(dyn_data.grid_v_c)), minimum(dyn_data.grid_p_full), minimum(dyn_data.grid_t_c), minimum(dyn_data.grid_tracers_c)
        @info "What a stupid code"
    end

    Time_Advance!(dyn_data)
  
    #@info "sec: ", integrator.time+1200, sum(abs.(grid_u_n)), sum(abs.(grid_v_n)), sum(abs.(grid_t_n)) , sum(abs.(grid_ps_n))
    #@info "max: ", maximum(abs.(grid_u_n)), maximum(abs.(grid_v_n)), maximum(abs.(grid_t_n)) , maximum(abs.(grid_ps_n))
    #@info "loc", grid_u_n[100,30,10],  grid_t_n[100,30,10], grid_u_n[1,32,1],  grid_t_n[1,32,1]
    
    #@assert(maximum(grid_u) <= 100.0 && maximum(grid_v) <= 100.0)


    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    
    return 

end 

function Get_Topography!(grid_geopots::Array{Float64, 3}, warm_start_file_name::String = "None", initial_day::Int64 = 5)
    # original_start = true
    # load_old_file  = false
    if warm_start_file_name == "None" # load warm start file
        grid_geopots .= 0.0
    end
    ### 2023/10/25
    # read_file     = load("/work/kaichiht/Colab/2023_research/annular_mode/300day_dry_run_all.dat")
    # grid_geopots .= read_file["grid_geopots_xyzt"][:,:,1,300]
    if warm_start_file_name != "None" # load warm start file
        read_file     = load(warm_start_file_name)
        grid_geopots .= read_file["grid_geopots_final"][:,:,:,1]
    end
    
    return
end 

function Spectral_Initialize_Fields!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, sea_level_ps_ref::Float64, init_t::Float64, 
                                     grid_geopots::Array{Float64,3}, dyn_data::Dyn_Data, Δt::Int64, warm_start_file_name::String = "None", initial_day::Int64 = 5)
    """
    Initializes the atmospheric state variables with a balanced background flow and specific perturbations.

    This function establishes the starting conditions for the simulation (t=0). It typically sets up a 
    zonally symmetric reference state (like a resting atmosphere or a jet) and superimposes small 
    perturbations to trigger baroclinic instability, which is essential for testing the dynamical 
    core's ability to develop cyclogenesis.

    The initialization sequence is:
    1.  **Thermodynamic Base State:**
        -   Sets a uniform initial temperature (`init_t`).
        -   Computes the surface log-pressure (`grid_lnps`) assuming hydrostatic balance with the 
            defined topography (`grid_geopots`) and a reference sea-level pressure.
    2.  **Spectral Perturbation:**
        -   Initializes vorticity and divergence to zero (resting air).
        -   Injects small-amplitude noise (`initial_perturbation`) into specific spherical harmonic 
            modes (e.g., wavenumber m=2, n=5) of the vorticity field. This breaks symmetry and 
            kickstarts wave growth.
    3.  **Consistent Wind Field:**
        -   Derives the grid-point wind field (`grid_u`, `grid_v`) from the perturbed spectral 
            vorticity and divergence to ensure consistency between spectral and grid representations.
    4.  **Spectral Transforms:**
        -   Performs forward and backward transforms on Temperature and Log-Pressure to ensure the 
            grid fields are spectrally truncated (filtered) correctly, removing any sub-grid scale 
            noise implied by the analytic setup.
    5.  **Pressure Variables:**
        -   Updates the 3D pressure thickness and interface variables (`grid_p_half`, `grid_Δp`) 
            consistent with the new surface pressure.

    Parameters
    ----------
    mesh (Spectral_Spherical_Mesh)
        Geometry and spectral transform infrastructure.

    atmo_data (Atmo_Data)
        Physical constants (Gas constant `rdgas`).

    vert_coord (Vert_Coordinate)
        Vertical grid structure for calculating pressure levels.

    sea_level_ps_ref (Float64)
        Reference pressure at mean sea level (e.g., 100,000 Pa).

    init_t (Float64)
        Isothermal temperature value for the initial state.

    grid_geopots (Array{Float64,3})
        Surface geopotential (topography).

    dyn_data (Dyn_Data)
        Data container where the initialized fields (u, v, T, ps, and their spectral counterparts) 
        will be stored.

    Returns
    -------
    None
        Operates in-place on `dyn_data`.

    Modified
    --------
    dyn_data.spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c
        Initialized spectral coefficients (current time step).

    dyn_data.grid_u_c, grid_v_c, grid_ps_c, grid_t_c
        Initialized grid-point fields (current time step).

    dyn_data.spe_X_p, grid_X_p
        Previous time step variables are initialized identically to the current step 
        (cold start condition).

    dyn_data.grid_p_half, grid_Δp, etc.
        Derived vertical pressure variables are updated.

    """

    if warm_start_file_name != "None" # load warm start file
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t           = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p   = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div            = dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd                                = mesh.nλ, mesh.nθ, mesh.nd
        
        ### By CJY2
        grid_t_n          = dyn_data.grid_t_n
        spe_tracers_c     = dyn_data.spe_tracers_c
        spe_tracers_p     = dyn_data.spe_tracers_p 

        grid_tracers_n    = dyn_data.grid_tracers_n
        grid_tracers_c    = dyn_data.grid_tracers_c
        grid_tracers_p    = dyn_data.grid_tracers_p 

        grid_u_n      = dyn_data.grid_u_n
        grid_v_n      = dyn_data.grid_v_n
        ########################################################
        # Tendency 
        grid_δu = dyn_data.grid_δu
        grid_δv = dyn_data.grid_δv

        grid_δtracers = dyn_data.grid_δtracers
        ########################################################
        @info warm_start_file_name # to make sure get the correct warmstart_PR.dat
        read_file     = load(warm_start_file_name)        
        grid_u[:,:,:]    .= read_file["grid_u_c_final"][:,:,:,1]
        grid_v[:,:,:]    .= read_file["grid_v_c_final"][:,:,:,1]  
        grid_t       .= read_file["grid_t_c_final"][:,:,:,1] 
        
        grid_lnps    .= log.(read_file["grid_ps_c_final"][:,:,1,1])
        grid_ps      .= read_file["grid_ps_c_final"][:,:,1,1]

        # grid_ps -> grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        ########################################################
        # By CJY
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical

        spe_t_c       .= read_file["spe_t_c_final"][:,:,:,1] 
        
        
        spe_vor_c[:,:,:] .= read_file["spe_vor_c_final"][:,:,:,1]
        spe_div_c[:,:,:] .= read_file["spe_div_c_final"][:,:,:,1]
        
        spe_lnps_c    .= (read_file["spe_lnps_c_final"][:,:,1,1])
        

        grid_vor .= read_file["grid_vor_final"][:,:,:,1] # Compute_Abs_Vor! need it 
        grid_div .= read_file["grid_div_final"][:,:,:,1] # Four_in_one! need it
        ########################################################
        spe_vor_p   .= read_file["spe_vor_p_final"][:,:,:,1]
        spe_div_p   .= read_file["spe_div_p_final"][:,:,:,1]
        spe_lnps_p  .= read_file["spe_lnps_p_final"][:,:,:,1]
        spe_t_p     .= read_file["spe_t_p_final"][:,:,:,1]

        grid_u_p    .= read_file["grid_u_p_final"][:,:,:,1]
        grid_v_p    .= read_file["grid_v_p_final"][:,:,:,1]
        grid_ps_p   .= read_file["grid_ps_p_final"][:,:,:,1]
        grid_t_p    .= read_file["grid_t_p_final"][:,:,:,1]
        ########################################################
        # Tracer initialization
        grid_tracers_n .= read_file["grid_tracers_n_final"][:,:,:,1] # large precipitation need next DO NOT REMOVE IT !!!
        grid_tracers_c .= read_file["grid_tracers_c_final"][:,:,:,1]
        grid_tracers_p .= read_file["grid_tracers_p_final"][:,:,:,1]
        
        # Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        # Trans_Grid_To_Spherical!(mesh, grid_tracers_p, spe_tracers_p)
        spe_tracers_c  .= read_file["spe_tracers_c_final"][:,:,:,1]
        spe_tracers_p  .= read_file["spe_tracers_p_final"][:,:,:,1]

        # grid_δu .= read_file["grid_δu_xyzt"][:,:,:,initial_day] # Rayleigh_Damping! has given it value 
        # grid_δv .= read_file["grid_δv_xyzt"][:,:,:,initial_day] # Rayleigh_Damping! has given it value

        ####################################################################
        # Correction_Init! would use these!!!
        # grid_t_n    .= read_file["grid_t_n_xyzt"][:,:,:,initial_day] 
        # grid_u_n   .= read_file["grid_u_n_xyzt"][:,:,:,initial_day]
        # grid_v_n   .= read_file["grid_v_n_xyzt"][:,:,:,initial_day]
        # grid_δtracers .= read_file["grid_δtracers_xyzt"][:,:,:,initial_day]

    end

    if warm_start_file_name == "None" # then use original start
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div =  dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
                
        ### By CJY2
        spe_tracers_c     = dyn_data.spe_tracers_c
        spe_tracers_p     = dyn_data.spe_tracers_p 
            
        grid_tracers_c    = dyn_data.grid_tracers_c
        grid_tracers_p    = dyn_data.grid_tracers_p 

        rdgas = atmo_data.rdgas
        grid_t    .= init_t

        # 𝚪 = 0.005
        # a = 6.371E6
        # b = 2
        # k = 3
        # p0 = 100000
        # Rd = 287
        # g  = 9.81
        # T0P = 240
        # T0E = 310
        # T0 = (T0E + T0P) * 0.5
        # H = Rd * T0/g
        
        # grid_z_full = zeros(((128,64,20)))
        # dry_run_file = load("test_final.dat")
        # grid_z_full .= dry_run_file["grid_z_full_xyzt"][:,:,:,5]
        
        # A = 1/𝚪 
        # B = (T0E - T0P) / (T0E + T0P) / T0P
        # C = (k+2)/2 * (T0E - T0P) / (T0E * T0P)
        
        # τ1 = zeros(((128,64,20)))
        # τ2 = zeros(((128,64,20)))
        
        # τ1 .= A * 𝚪 / T0 .* exp.(𝚪/T0 .* grid_z_full) .+ B .* (1 .- 2 .* (grid_z_full./(b*H)).^2) .* exp.(-1 .* (grid_z_full./(b*H)).^2) 
        
        # τ2 .= C .* (1 .- 2 .* (grid_z_full./(b*H)).^2) .* exp.(-1 .* (grid_z_full./(b*H)).^2) 
        
        # θc2  = LinRange(-90,90,64)
        # θc  = deg2rad.(θc2)
        # for j in 1:64
        #     grid_t[:,j,:] .= (τ1[:,j,:] .- τ2[:,j,:] .* ((cos(θc[j]))^k - (k/(k+2)) .* (cos(θc[j]))^(k+2))).^-1
        # end

        
        # dΦ/dlnp = -RT    Δp = -ΔΦ/RT
        grid_lnps[:,:,1] .= log(sea_level_ps_ref) .- grid_geopots[:,:,1] ./ (rdgas * init_t) 
        grid_ps   .= exp.(grid_lnps)
        
        
        spe_div_c .= 0.0
        spe_vor_c .= 0.0
      
        # # initial perturbation
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
        
        initial_perturbation = 1.0e-7/sqrt(2.0)
        # initial vorticity perturbation used in benchmark code
        # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
        for k = nd-2:nd
          spe_vor_c[2,5,k] = initial_perturbation
          spe_vor_c[6,9,k] = initial_perturbation
          spe_vor_c[2,4,k] = initial_perturbation  
          spe_vor_c[6,8,k] = initial_perturbation
        end
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
      
        # initial spectral fields (and spectrally-filtered) grid fields
        
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
      
        Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
        Trans_Spherical_To_Grid!(mesh, spe_lnps_c,  grid_lnps)
        grid_ps .= exp.(grid_lnps)
        
      
        Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
        Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
        Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)
        
        #update pressure variables for hs forcing
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        
        
        spe_vor_p  .= spe_vor_c
        spe_div_p  .= spe_div_c
        spe_lnps_p .= spe_lnps_c
        spe_t_p    .= spe_t_c
      
      
        grid_u_p   .= grid_u
        grid_v_p   .= grid_v
        grid_ps_p  .= grid_ps
        grid_t_p   .= grid_t

        # Tracer initialization
        initial_RH      = 0.8
        Lv              = 2.5*10^6.
        Rv              = 461.
        qv0             = 0.018
        θc              = mesh.θc # lat
        phi_hw          = 2 * pi / 9 * deg2rad(40)
        p_hw            = 30000.
        phi             = LinRange(-90,90,64)
        p0              = 100000.
        for k in 1:20
            for j in 1:64
               for i in 1:128
                   # grid_tracers_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((deg2rad(phi[j]))/phi_hw)^4) 
                   grid_tracers_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((θc[j])/phi_hw)^4) 
                    
               end            
            end
        end
        grid_tracers_c[:,:,1] .= 0.
        
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)

        grid_tracers_p .= grid_tracers_c
        spe_tracers_p  .= spe_tracers_c
    end
      

     
end 


function Spectral_Dynamics_Physics!(semi_implicit::Semi_Implicit_Solver, atmo_data::Atmo_Data, mesh::Spectral_Spherical_Mesh, dyn_data::Dyn_Data, Δt::Int64, physics_params::Dict{String, Float64}, L::Float64)
    """
    Computes the physical forcing tendencies for the atmospheric state variables.

    This function acts as the interface between the dynamical core and the physical parameterizations.
    In the current configuration, it implements the **Held-Suarez (1994)** benchmark forcing, which is
    standard for testing the climatology of dynamical cores without the complexity of moist physics or
    radiation.

    The function performs the following operations:
    1.  **Initialization:** Resets the surface pressure tendency (`grid_δps`) to zero, as the
        simplified physics does not directly add mass sources/sinks (unlike precipitation in a full
        model).
    2.  **Held-Suarez Forcing:** Calls `HS_Forcing!` to calculate:
        -   **Thermal Relaxation:** Relaxes the temperature field (`grid_t`) towards a prescribed
            zonally symmetric radiative equilibrium state (`grid_t_eq`) using Newtonian cooling.
        -   **Rayleigh Friction:** Applies a linear damping to the horizontal winds (`grid_u`,
            `grid_v`) near the surface to simulate boundary layer friction.
    3.  **Tendency Update:** The calculated physical tendencies are added directly to the accumulation
        arrays (`grid_δt`, `grid_δu`, `grid_δv`).

    Parameters
    ----------
    atmo_data (Atmo_Data)
        Structure containing physical constants and model configuration.

    mesh (Spectral_Spherical_Mesh)
        Geometry information (specifically `sinθ` for latitude-dependent forcing).

    dyn_data (Dyn_Data)
        Container holding the prognostic variables (predictor state `_p`) and tendency arrays.

    Δt (Int64)
        The physics time step in seconds.

    physics_params (Dict{String, Float64})
        A dictionary of coefficients controlling the forcing (e.g., relaxation time scales `ka`, `ks`,
        or friction coefficient `kf`).

    Returns
    -------
    None
        The function operates in-place on the tendency arrays within `dyn_data`.

    Modified
    --------
    dyn_data.grid_δu, dyn_data.grid_δv
        Decremented by Rayleigh friction (mostly in the lower atmosphere).

    dyn_data.grid_δt
        Incremented/Decremented by Newtonian cooling (thermal relaxation).

    dyn_data.grid_δps
        Reset to 0.0 (explicitly handling the mass budget for physics).
    """
    grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    
    # spectral quantities
    spe_t_c       = dyn_data.spe_t_c
    spe_tracers_c = dyn_data.spe_tracers_c

    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t_c, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n

    # related quanties
    grid_p_half, grid_p_full = dyn_data.grid_p_half, dyn_data.grid_p_full
    grid_t_eq = dyn_data.grid_t_eq

    # By CJY2
    grid_tracers_p, grid_tracers_c, grid_tracers_n    = dyn_data.grid_tracers_p, dyn_data.grid_tracers_c, dyn_data.grid_tracers_n
    grid_tracers_diff = dyn_data.grid_tracers_diff
    grid_δtracers     = dyn_data.grid_δtracers 

    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)
    factor1           = dyn_data.factor1 
    factor2           = dyn_data.factor2 
    factor3           = dyn_data.factor3  
    # factor4 = dyn_data.factor4  
    K_E               = dyn_data.K_E
    ###
    
    grid_δps .= 0.0
    grid_δtracers  .= 0.

    ### Cal V_c and za
    V_c, za, rho = Calculate_V_c_za_rho(atmo_data, dyn_data, grid_p_half, grid_p_full, grid_ps, grid_u, grid_v, grid_t_c, grid_tracers_c)
    
    """
    ## large-scale precipitation
    """
    do_large_scale_precipitation = true
    do_Sensible_heat_fluxes      = true
    do_Surface_evaporation       = true
    do_Implicit_PBL_Scheme       = true
    
    if do_large_scale_precipitation == true
        lscale_cond!(
            semi_implicit, grid_p_full,
            grid_t_c, grid_tracers_c,
            grid_δt, grid_δtracers, factor3, 
            L
        )

        grid_tracers_c[grid_tracers_c .< 0]   .= 0     
    
        
        grid_tracers_c .= grid_tracers_c .- grid_δtracers .* (2*Δt)
        grid_t_c       .= grid_t_c       .+ grid_δt       .* (2*Δt)
    
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
        
        Trans_Grid_To_Spherical!(mesh, grid_t_c, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t_c)
        
        grid_δtracers .= 0.
        grid_δt       .= 0.
    end

    # Calculate grid_δt(.+=) and grid_t(.=)
    if do_Sensible_heat_fluxes == true
        Sensible_heating!(
            mesh, atmo_data, 
            grid_t_c, 
            V_c, za, Δt
        )
        Trans_Grid_To_Spherical!(mesh, grid_t_c, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t_c)
    end

    if do_Surface_evaporation == true
        # Calculate grid_δtracers(.+=) and grid_tracers_c(.=)  (Latent_heat_flux! == Surface_evaporation!)
        Surface_evaporation!(
            mesh, atmo_data,
            grid_ps,
            grid_tracers_c,
            V_c, za, Δt, factor1
        )
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
    end

    # Calculate {grid_δtracers(.+=) and grid_tracers_c(.=)} and {grid_δt(.+=) and grid_t(.=)}
    if do_Implicit_PBL_Scheme == true
        Implicit_PBL_Scheme!(
            atmo_data,
            grid_t_c, grid_tracers_c,
            grid_δt, grid_δtracers,
            grid_p_full, grid_p_half,
            K_E, V_c, za, rho,
            Δt, factor2
        )
        Trans_Grid_To_Spherical!(mesh, grid_t_c, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t_c)
        
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
    end
    
    # Note: grid_δt was zeroed before this step
    ######################################################################################################
    # summer 2025
    # edit: pass 2 more args to HS_Forcing!
    # "grid_tracers_p" for LRF calc
    # "dyn_data.grid_δt_HS" for individual forcing output
    HS_Forcing!(
        atmo_data, Δt, mesh.sinθ, 
        grid_u_p, grid_v_p, 
        grid_p_half, grid_p_full, 
        grid_t_p, grid_tracers_p,
        grid_δu, grid_δv,
        grid_t_eq, grid_δt, dyn_data.grid_δt_HS, 
        physics_params
    )

end


function Atmosphere_Update!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, semi_implicit::Semi_Implicit_Solver, 
                            dyn_data::Dyn_Data, physcis_params::Dict{String, Float64}, L::Float64)

    Δt = Get_Δt(semi_implicit.integrator)
    Spectral_Dynamics_Physics!(semi_implicit, atmo_data, mesh,  dyn_data, Δt, physcis_params, L) # HS forcing
    Spectral_Dynamics!(mesh,  vert_coord , atmo_data, dyn_data, semi_implicit) # dynamics 

    grid_ps , grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_ps_c,  dyn_data.grid_Δp, dyn_data.grid_p_half, 
                                                                                dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full 
    
    grid_t = dyn_data.grid_t_c
    grid_geopots, grid_z_full, grid_z_half = dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

    ### 1201
    grid_tracers_c = dyn_data.grid_tracers_c
        
    Compute_Pressures_And_Heights!(atmo_data, vert_coord,     
    grid_ps, grid_geopots, grid_t, 
    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full, grid_z_full, grid_z_half, grid_tracers_c)

    return
end 



# function Spectral_Dynamics_Main()
#   # the decay of a sinusoidal disturbance to a zonally symmetric flow 
#   # that resembles that found in the upper troposphere in Northern winter.
#   name = "Spectral_Dynamics"
#   #num_fourier, nθ, nd = 63, 96, 20
#   num_fourier, nθ, nd = 42, 64, 20
#   #num_fourier, nθ, nd = 21, 32, 20
#   num_spherical = num_fourier + 1
#   nλ = 2nθ
  
#   radius = 6371000.0
#   omega = 7.292e-5
#   sea_level_ps_ref = 1.0e5
#   init_t = 264.0
  
#   # Initialize mesh
#   mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
#   θc, λc = mesh.θc,  mesh.λc
#   cosθ, sinθ = mesh.cosθ, mesh.sinθ
  
#   vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
#   # Initialize atmo_data
#   do_mass_correction = true
#   do_energy_correction = true
#   do_water_correction = false
  
#   use_virtual_temperature = false
#   atmo_data = Atmo_Data(name, nλ, nθ, nd, do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature, sinθ, radius,  omega)
  
#   # Initialize integrator
#   damping_order = 4
#   damping_coef = 1.15741e-4
#   robert_coef  = 0.04 
  
#   implicit_coef = 0.5
#   day_to_sec = 86400
#   start_time = 0
#   end_time = 2*day_to_sec  #
#   Δt = 1200
#   init_step = true
  
#   integrator = Filtered_Leapfrog(robert_coef, 
#   damping_order, damping_coef, mesh.laplacian_eig,
#   implicit_coef, Δt, init_step, start_time, end_time)
  
#   ps_ref = sea_level_ps_ref
#   t_ref = fill(300.0, nd)
#   wave_numbers = mesh.wave_numbers
#   semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data,
#   integrator, ps_ref, t_ref, wave_numbers)
  
#   # Initialize data
#   dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)
  
  
#   NT = Int64(end_time / Δt)
  
#   Get_Topography!(dyn_data.grid_geopots)
  
#   Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t,
#   dyn_data.grid_geopots, dyn_data)
  

#   Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#   Update_Init_Step!(semi_implicit)
#   integrator.time += Δt 
#   for i = 2:NT

#     Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#     integrator.time += Δt
#     @info integrator.time

#   end
  
# end


# #Spectral_Dynamics_Main()

