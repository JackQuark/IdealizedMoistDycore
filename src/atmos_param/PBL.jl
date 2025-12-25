# Quark 12/25,2025

function Calculate_V_c_za_rho(
    atmo_data::Atmo_Data, dyn_data::Dyn_Data,
    grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_ps::Array{Float64, 3},
    grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3},
    grid_t::Array{Float64, 3}, grid_q::Array{Float64, 3}
)
    Lv  = atmo_data.Lv
    Rv  = atmo_data.rvgas
    Rd  = atmo_data.rdgas
    cp  = atmo_data.cp_air
    grav = atmo_data.grav

    # horizontal wind speed
    V_c = sqrt.(grid_u .^ 2 .+ grid_v .^ 2)
    
    # height of the lowest model level
    tv = grid_t[:,:,20] .* (1. .+ 0.608 .* grid_q[:,:,20])
    za = Rd .* tv[:,:] ./grav .* (log.(grid_ps[:,:,1] ./ ((grid_p_full[:,:,20] .+ grid_p_half[:,:,21]) ./ 2.) )) ./2

    # density
    rho = grid_p_full[:,:,:] ./ Rd ./ (grid_t[:,:,:].* (1. .+ 0.608 .* grid_q[:,:,20]))
    
    return V_c, za, rho
end

"""heating directly operates on `grid_t`"""
function Sensible_heating!(
    mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_t::Array{Float64, 3},
    V_c::Array{Float64, 3}, za::Array{Float64, 2}, Δt::Int64,
    C_H::Float64=0.0044
)
    θc  = mesh.θc
    Tsurf = zeros((128,64))
    for i in 1:64
         Tsurf[:,i] .= 29. .* exp.(-(θc[i] .^2. ./ (2 * (26. * pi / 180.)^2.))) .+ 271.
    end

    grid_t[:,:,20]  .= (grid_t[:,:,20] .+ C_H .* V_c[:,:,20] .* Tsurf[:,:,1] .* Δt ./ za[:,:]) ./ 
                       (1. .+ C_H .* V_c[:,:,20] .* Δt ./ za[:,:])

end

"""evaporation directly operates on `grid_q`"""
function Surface_evaporation!(
    mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_ps::Array{Float64, 3},
    grid_q::Array{Float64,3},
    V_c::Array{Float64, 3}, za::Array{Float64, 2},
    Δt::Int64, factor1::Array{Float64, 3},
    C_E::Float64=0.0044
)
    Lv  = atmo_data.Lv
    Rv  = atmo_data.rvgas
    Rd  = atmo_data.rdgas
    θc = mesh.θc
    Tsurf = zeros((128,64))
    for i in 1:64
        Tsurf[:,i] .= 29. .* exp.(-(θc[i] .^2. ./ (2 * (26. * pi / 180.)^2.))) .+ 271.
    end
    
    grid_esat_ps = 611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ Tsurf[:,:]))
    grid_qsat_ps = (0.622 .* grid_esat_ps) ./ (grid_ps[:,:,1] .- 0.378 .* grid_esat_ps)

    grid_q[:,:,20]      .= ((grid_q[:,:,20] .+ C_E .* V_c[:,:,20] .* max.(grid_q[:,:,20], grid_qsat_ps[:,:,1]) .* Δt ./ za[:,:]) ./ (1. .+ C_E .* V_c[:,:,20]  .* Δt ./ za[:,:]))

    factor1[:,:,20]     .= grid_q[:,:,20] ./(2. .* Δt) # TODO: fix bad naming
    
end 

function Implicit_PBL_Scheme!(
    atmo_data::Atmo_Data,
    grid_t::Array{Float64, 3}, grid_q::Array{Float64, 3},
    grid_δt::Array{Float64, 3}, grid_δq::Array{Float64, 3},
    grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3},
    K_E::Array{Float64, 3}, V_c::Array{Float64, 3}, za::Array{Float64, 2}, rho::Array{Float64, 3},
    Δt::Int64, factor2::Array{Float64, 3},
    C_D::Float64=0.0044
)
    """
    Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
    
    Step1. Calculate K_E = C_E * V_c * za
    Step2. Calculate A, B, C --> E, F
    Step3. Calculate new grid_q and grid_t
    """
    Lv  = atmo_data.Lv
    Rv  = atmo_data.rvgas
    Rd  = atmo_data.rdgas
    cp  = atmo_data.cp_air
    
    V_a = V_c[:,:,20]
    
    grav = atmo_data.grav
    
    ### 12/28 upgrade output_manager
    K_E[:,:,17:20] .= C_D .* V_a .* za[:,:]
    K_E[:,:, 1:16] .= C_D .* V_a .* za[:,:] .* exp.(-((85000. .- grid_p_half[:,:,1:16]) ./ 10000.).^2)
    ### cal PBL Scheme
    rpdel  = zeros(((128,64,20))) ### = 1 / (p^n_{+} - p^n_{-}) , which p^_{-} mean upper layer
    for i in 1:20
        rpdel[:,:,i] .= 1. ./ (grid_p_half[:,:,i+1] .- grid_p_half[:,:,i])
    end

    CA     = zeros(((128,64,20)))
    CC     = zeros(((128,64,20)))
    CE     = zeros(((128,64,20+1)))
    CF     = zeros(((128,64,20+1)))
    CFt    = zeros(((128,64,20+1)))

    # replace for loop by array operation >>>
    ks  = 1:19
    k1s = 2:20
    Δp_full = grid_p_full[:, :, k1s] .- grid_p_full[:, :, ks]
    CA[:, :, ks]  .= rpdel[:, :, ks]  .* Δt .* grav^2 .* K_E[:, :, k1s] .* rho[:, :, k1s].^2 ./ Δp_full
    CC[:, :, k1s] .= rpdel[:, :, k1s] .* Δt .* grav^2 .* K_E[:, :, k1s] .* rho[:, :, k1s].^2 ./ Δp_full
    # <<< modified by Quark 12/25,2025

    CA[:,:,20]   .= 0.
    CC[:,:, 1]   .= 0.
    CE[:,:, 1]   .= 0.
    CE[:,:,21]   .= 0.
    CF[:,:,21]   .= 0.
    CFt[:,:,21]  .= 0.
    
    p0 = 100000.
    for k in 20:-1:1
        CE[:,:,k]      .= CC[:,:,k] ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1])

        CF[:,:,k]      .= ((grid_q[:,:,k] .+ CA[:,:,k] .* CF[:,:,k+1]) ./
                           (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1]))

        CFt[:,:,k]     .= (((p0./grid_p_full[:,:,k]).^(Rd/cp).*grid_t[:,:,k] .+ CA[:,:,k] .* CFt[:,:,k+1]) ./
                           (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1]))
    end

    # first calculate the updates at the top model level
    # grid_δq[:,:,1] .+= (CF[:,:,1] .- grid_q[:,:,1]) ./ (2. .* Δt)
    ### WARNING factor1 just factor, so it did  ./ ./ (2. .* Δt). 
    ### So did factor2
    factor2[:,:,1]     .= (CF[:,:,1] .- grid_q[:,:,1]) ./ (2. .* Δt) # because CE at top = 0 # TODO: fix bad naming
    grid_q[:,:,1]      .= CF[:,:,1] 
    grid_t[:,:,1]      .= (CFt[:,:,1] .* (grid_p_full[:,:,1]./p0).^(Rd/cp))
    
    # Loop over the remaining level
    for k in 2:20
        factor2[:,:,k] .= (CE[:,:,k] .* grid_q[:,:,k-1] .+ CF[:,:,k] .- grid_q[:,:,k]) ./ (2. .* Δt)
        grid_q[:,:,k]  .=  CE[:,:,k] .* grid_q[:,:,k-1] .+ CF[:,:,k]
        grid_t[:,:,k]  .= ((CE[:,:,k] .* grid_t[:,:,k-1] .* (p0./grid_p_full[:,:,k-1]).^(Rd/cp) .+ CFt[:,:,k]) .* (grid_p_full[:,:,k]./p0).^(Rd/cp))
    end

end