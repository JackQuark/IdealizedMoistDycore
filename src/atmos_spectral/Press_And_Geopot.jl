export Compute_Pressures_And_Heights!, Half_Level_Pressures!, Pressure_Variables!, Compute_Geopotential!

function Compute_Pressures_And_Heights!(atmo_data::Atmo_Data, vert_coord::Vert_Coordinate,     
                                        grid_ps::Array{Float64,3}, grid_geopots::Array{Float64,3}, grid_t::Array{Float64,3}, 
                                        grid_p_half::Array{Float64,3},  grid_Δp::Array{Float64,3},
                                        grid_lnp_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, grid_lnp_full::Array{Float64,3},
                                        grid_z_full::Array{Float64,3}, grid_z_half::Array{Float64,3}, grid_tracers_c::Array{Float64,3})
    """
    Updates the diagnostic pressure and geometric height variables for the entire atmospheric column,
    ensuring hydrostatic balance consistent with the provided vertical coordinate system.

    This function performs the following computational steps:
    1. Reconstructs 3D pressure fields on both layer centers (full) and interfaces (half) using
    surface pressure and vertical grid coefficients.
    2. Computes the logarithmic pressure thickness (Δlnp) and layer thickness (Δp) required for
    vertical integration.
    3. Integrates the hydrostatic equation to solve for geopotential (Φ) at all levels, utilizing the
    provided temperature field.
    4. Normalizes geopotential by the gravitational constant (g) to update the geometric height (z)
    arrays for both full and half levels.

    Parameters
    ----------
    atmo_data (Atmo_Data)
        Structure containing physical constants for the simulation, specifically the gravitational
        acceleration (grav).

    vert_coord (Vert_Coordinate)
        Vertical discretization definitions, including hybrid coefficients (A_k, B_k) used to map
        surface pressure to 3D pressure levels.

    grid_ps (Array{Float64,3})
        Surface pressure field (p_s) used as the boundary condition for vertical integration.

    grid_geopots (Array{Float64,3})
        Surface geopotential field (Φ_s), representing surface topography.

    grid_t (Array{Float64,3})
        Atmospheric temperature (T) on full vertical levels, used in the hydrostatic integration.
    
    Returns
    -------
    None
        The function operates in-place, modifying the provided pressure and height arrays directly.
    
    Modified
    --------
    grid_p_half (Array{Float64,3})
        Overwritten with pressure values at layer interfaces (p_{k+1/2}).

    grid_Δp (Array{Float64,3})
        Overwritten with pressure thickness of each layer (Δp_k).

    grid_lnp_half (Array{Float64,3})
        Overwritten with the natural logarithm of interface pressures (ln p_{k+1/2}).

    grid_p_full (Array{Float64,3})
        Overwritten with pressure values at layer centers (p_k).

    grid_lnp_full (Array{Float64,3})
        Overwritten with the natural logarithm of layer center pressures (ln p_k).

    grid_z_full (Array{Float64,3})
        Overwritten with geometric height at layer centers (z_k).

    grid_z_half (Array{Float64,3})
        Overwritten with geometric height at layer interfaces (z_{k+1/2}).

    """
    grav = atmo_data.grav

    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
                        grid_lnp_half, grid_p_full, grid_lnp_full)

    Compute_Geopotential!(vert_coord,atmo_data, grid_lnp_half, grid_lnp_full,  grid_t, grid_geopots, grid_z_full, grid_z_half, grid_tracers_c)

    grid_z_full ./= grav
    grid_z_half ./= grav

end 



function Half_Level_Pressures!(vert_coord::Vert_Coordinate, grid_ps::Array{Float64,3}, grid_p_half::Array{Float64,3}) 
    """
    Computes the hydrostatic pressure at vertical layer interfaces (half-levels)
    based on the hybrid sigma-pressure coordinate definition.

    This function performs the following computational steps:
    1. Accesses the hybrid vertical coordinate coefficients (A_k and B_k) and the number of vertical levels (nd) from the coordinate structure.
    2. Iterates through each vertical interface index (k from 1 to nd+1).
    3. Calculates the pressure at each interface using the linear relationship: p_{k+1/2} = A_k + B_k * p_s, where p_s is the surface pressure.

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        Structure defining the vertical discretization, containing the hybrid coefficients 'ak' (static pressure component) and 'bk' (sigma component).

    grid_ps (Array{Float64,3})
        Surface pressure field (p_s), serving as the dynamic lower boundary condition. Dimensions are typically (Longitude, Latitude, 1).

    grid_p_half (Array{Float64,3})
        Buffer to be updated with the computed pressure values at layer interfaces (p_{k+1/2}). Dimensions are typically (Longitude, Latitude, Levels+1).

    Returns
    -------
    None
        The function operates in-place, populating the 'grid_p_half' array.
    
    Modified
    --------
    grid_p_half (Array{Float64,3})
        Overwritten with the computed pressure values at layer interfaces (p_{k+1/2}). Dimensions are typically (Longitude, Latitude, Levels+1).
    
    """
    nd = vert_coord.nd
    bk = vert_coord.bk
    ak = vert_coord.ak
  
    # pk = ak * pref
    for k=1:nd+1
        grid_p_half[:,:,k] .= ak[k] .+ bk[k]*grid_ps[:,:,1]
    end 
end 


function Pressure_Variables!(vert_coord::Vert_Coordinate, grid_ps::Array{Float64,3}, grid_p_half::Array{Float64,3}, grid_Δp::Array{Float64,3},
                             grid_lnp_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, grid_lnp_full::Array{Float64,3})
    """
    Computes the full suite of diagnostic pressure variables required for the dynamical core,
    deriving full-level and interface values from surface pressure.

    This function implements the energy-conserving vertical discretization scheme described by
    Simmons and Burridge (1981). It performs the following operations:
    1. Validates input dimensions (surface pressure must be 2D).
    2. Calls `Half_Level_Pressures!` to compute interface pressures (p_{k+1/2}).
    3. Calculates the pressure thickness (Δp) for every layer.
    4. Computes the natural logarithm of pressures at interfaces.
    5. Derives the logarithm of pressure at layer centers (ln p_k) using the Simmons-Burridge
    formulation, which ensures analytical consistency for angular momentum and energy conservation.
    Special handling is applied for the top-most layer if the model top is at zero pressure.
    6. Recovers the full-level pressure (p_k) by exponentiating the derived log values.

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        Configuration for the vertical grid, including `vert_difference_option` (must be
        "simmons_and_burridge") and the `zero_top` flag indicating if the top interface pressure is
        strictly zero.

    grid_ps (Array{Float64,3})
        Surface pressure field (p_s). The third dimension must be singleton (size 1).

    grid_p_half (Array{Float64,3})
        Buffer for pressure values at layer interfaces (p_{k+1/2}).

    grid_Δp (Array{Float64,3})
        Buffer for the pressure thickness of each layer (Δp_k).

    grid_lnp_half (Array{Float64,3})
        Buffer for the natural logarithm of interface pressures (ln p_{k+1/2}).

    grid_p_full (Array{Float64,3})
        Buffer for pressure values at layer centers (p_k).

    grid_lnp_full (Array{Float64,3})
        Buffer for the natural logarithm of layer center pressures (ln p_k).

    Returns
    -------
    None
        The function operates in-place, populating the provided pressure arrays via the modified
        arguments.

    Modified
    --------
    grid_p_half
        Overwritten with computed interface pressures.

    grid_Δp
        Overwritten with computed layer pressure thicknesses.

    grid_lnp_half
        Overwritten with computed log interface pressures.

    grid_p_full
        Overwritten with computed layer center pressures.

    grid_lnp_full
        Overwritten with computed log layer center pressures.
        
    """
    @assert(size(grid_ps)[3] == 1)
    
    Half_Level_Pressures!(vert_coord, grid_ps, grid_p_half)
    nd = vert_coord.nd
    grid_Δp .= grid_p_half[:,:,2:nd+1] - grid_p_half[:,:,1:nd]
    zero_top = vert_coord.zero_top
    
    if (vert_coord.vert_difference_option == "simmons_and_burridge") 
    
        k_top = (zero_top ? 2 : 1) 
        
        grid_lnp_half[:,:,k_top:nd+1] .= log.(grid_p_half[:,:,k_top:nd+1])
        
        #lnp_{k} = (p_{k+1/2}lnp_{k+1/2} - p_{k-1/2}lnp_{k-1/2})/Δp_k - 1
        #        = [(p_{k+1/2}-p_{k-1/2})lnp_{k+1/2} + p_{k-1/2}(lnp_{k+1/2} - lnp_{k-1/2})]/Δp_k - 1
        #        = lnp_{k+1/2} + [p_{k-1/2}(lnp_{k+1/2} - lnp_{k-1/2})]/Δp_k - 1
        grid_lnp_full[:,:,k_top:nd] .= grid_lnp_half[:,:,k_top+1:nd+1] .+ grid_p_half[:,:,k_top:nd].*(grid_lnp_half[:,:,k_top+1:nd+1] - grid_lnp_half[:,:,k_top:nd])./grid_Δp[:,:,k_top:nd] .- 1.0 
        
        if (zero_top) 
            grid_lnp_half[:,:,1] .= 0.0
            grid_lnp_full[:,:,1] .= grid_lnp_half[:,:,2] .- 1.0 
        end
    
    else
        error("vert_difference_option ",vert_coord.vert_difference_option, " is not a valid value for option")
    
    end
    grid_p_full .= exp.(grid_lnp_full)
end



function Compute_Geopotential!(vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
                               grid_lnp_half::Array{Float64, 3}, grid_lnp_full::Array{Float64, 3},  
                               grid_t::Array{Float64, 3}, 
                               grid_geopots::Array{Float64, 3}, grid_geopot_full::Array{Float64, 3}, grid_geopot_half::Array{Float64, 3},
                               grid_tracers_c::Array{Float64, 3})
    """
    Integrates the hydrostatic relation vertically to compute the geopotential field (Φ) at both layer
    interfaces and layer centers.

    This function performs a bottom-up integration starting from the surface topography:
    1. Initializes the bottom interface (k = nd+1) with the surface geopotential (grid_geopots).
    2. Integrates upward to compute geopotential at all interfaces (half-levels) using the provided
    temperature and logarithmic pressure gradients.
    3. Interpolates/Computes geopotential at layer centers (full-levels) consistent with the
    half-level values and the layer's temperature.

    Current implementation notes:
    - Uses the dry air gas constant (rdgas).
    - Virtual temperature correction is currently disabled (virtual_t = grid_t), assuming a dry
    atmosphere.

    Parameters
    ----------
    vert_coord (Vert_Coordinate)
        Vertical grid configuration, including the number of levels (nd) and the `zero_top` flag for
        handling the top-most boundary condition.

    atmo_data (Atmo_Data)
        Physical constants structure providing the gas constant for dry air (rdgas) and flags for
        thermodynamics.

    grid_lnp_half (Array{Float64, 3})
        Natural logarithm of pressure at layer interfaces (ln p_{k+1/2}).

    grid_lnp_full (Array{Float64, 3})
        Natural logarithm of pressure at layer centers (ln p_k).

    grid_t (Array{Float64, 3})
        Atmospheric temperature (T) on full vertical levels.

    grid_geopots (Array{Float64, 3})
        Surface geopotential height (Φ_s), serving as the lower boundary condition for the integration.

    grid_geopot_full (Array{Float64, 3})
        Buffer for geopotential values at layer centers (Φ_k).

    grid_geopot_half (Array{Float64, 3})
        Buffer for geopotential values at layer interfaces (Φ_{k+1/2}).

    Returns
    -------
    None
        The function operates in-place, populating the provided geopotential arrays via the modified
        arguments.

    Modified
    --------
    grid_geopot_full
        Overwritten with computed geopotential at layer centers.

    grid_geopot_half
        Overwritten with computed geopotential at layer interfaces.

    """
    use_virtual_temperature = atmo_data.use_virtual_temperature
    rvgas, rdgas = atmo_data.rvgas, atmo_data.rdgas
    zero_top = vert_coord.zero_top
    nd = vert_coord.nd


    grid_geopot_half[:,:,nd+1] .= grid_geopots[:,:,1]
    
    if zero_top  #todo (pk(1).eq.0.0) then
        k_top = 2
        grid_geopot_half[:,:,1] .= 0.0
    else
        k_top = 1
    end
    
    if (use_virtual_temperature) 
        virtual_t = grid_t .* (1. .+ (rvgas/rdgas - 1.).*grid_tracers_c)
    else
        virtual_t = grid_t
    end


    # virtual_t = grid_t
    
    
    for k=nd:-1:k_top
        #Φ_{k-1/2} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k-1})
        grid_geopot_half[:,:,k] .= grid_geopot_half[:,:,k+1] .+ rdgas*virtual_t[:,:,k] .* (grid_lnp_half[:,:,k+1] - grid_lnp_half[:,:,k])
    end
    
    for k=1:nd
        #Φ_{k} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k})
        grid_geopot_full[:,:,k] .= grid_geopot_half[:,:,k+1] .+ rdgas*virtual_t[:,:,k] .* (grid_lnp_half[:,:,k+1] - grid_lnp_full[:,:,k])
    end
  
end 
