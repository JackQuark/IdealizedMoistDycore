# edit: arg "grid_qv" added
function HS_Forcing!(
    atmo_data::Atmo_Data, Δt::Int64, sinθ::Array{Float64, 1}, 
    grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3}, 
    grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, 
    grid_t::Array{Float64, 3}, grid_qv::Array{Float64, 3}, 
    grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, 
    grid_t_eq::Array{Float64, 3}, grid_δt::Array{Float64, 3}, grid_δt_HS::Dict{String, Array{Float64, 3}},
    physcis_params::Dict{String, Float64}
    )

    σ_b = physcis_params["σ_b"]   
    k_f = physcis_params["k_f"]   # day^{-1}
    k_a = physcis_params["k_a"]   # day^{-1}
    k_s = physcis_params["k_s"]   # day^{-1}
    ΔT_y = physcis_params["ΔT_y"] # K
    Δθ_z = physcis_params["Δθ_z"] # K

    # rayleigh damping of wind components near the surface
    Rayleigh_Damping!(atmo_data, grid_p_half, grid_p_full, grid_u, grid_v, grid_δu, grid_δv, σ_b, k_f)

    #todo 
    do_conserve_energy = true
    if (do_conserve_energy) 
        grid_δt .= -((grid_u + 0.5*grid_δu*Δt).*grid_δu + (grid_v + 0.5*grid_δv*Δt).*grid_δv)/atmo_data.cp_air
    end
    # !!! forcing below sould be added to grid_δt_HS.
    # !!! Use self defined key name for each forcing,
    # !!! and add the key to the "keys" below.
    # !!! pls overwrite the existing Array in grid_δt_HS.
    keys = [
        "HS_Newtonian",
        # "LRF"
        ]
    for k in keys
        if !haskey(grid_δt_HS, k)
            grid_δt_HS[k] = zeros(Float64, size(grid_t))
        end
    end
    
    #thermal forcing for held & suarez (1994) benchmark calculation
    Newtonian_Damping!(atmo_data, sinθ, grid_p_half, grid_p_full, grid_t, grid_t_eq, grid_δt_HS, σ_b, k_a, k_s, ΔT_y, Δθ_z)
    
    # todo: LRF, t|qv perturbation (summer2025)
    # Radiative_Forcing!(
    #     atmo_data, Δt, grid_p_full, grid_t, grid_qv, grid_δt_HS,
    #     if_t_LRF=true, if_q_LRF=true
    # )

    for k in keys
        grid_δt .+= grid_δt_HS[k]
    end
end 


function Newtonian_Damping!(atmo_data::Atmo_Data, sinθ::Array{Float64,1}, 
    grid_p_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, 
    grid_t::Array{Float64,3}, grid_t_eq::Array{Float64,3}, grid_δt_HS::Dict{String, Array{Float64, 3}},
    σ_b::Float64, k_a::Float64, k_s::Float64, ΔT_y::Float64, Δθ_z::Float64)
    
    #routine to compute thermal forcing for held & suarez (1994)
    nλ, nθ, nd = size(grid_t)
    
    day_to_sec = 86400
    t_zero, t_strat =  315.0, 200.0
    k_a, k_s = k_a/day_to_sec,  k_s/day_to_sec#s^-1
    
    p_ref = 1.0e5

    grid_ps = grid_p_half[:,:,nd+1]
    
    sinθ_2d = repeat(sinθ', nλ, 1)
    sinθ2_2d = sinθ_2d .* sinθ_2d
    cosθ2_2d = 1.0  .- sinθ2_2d
    cosθ4_2d = cosθ2_2d .* cosθ2_2d
    
    kappa = atmo_data.kappa
    grid_p_norm  = zeros(Float64, nλ, nθ)
    σ = zeros(Float64, nλ, nθ)
    k_t = zeros(Float64, nλ, nθ)
    # compute equilibrium temperature (grid_t_eq)
    
    for  k = 1:nd
        grid_p_norm .= grid_p_full[:,:,k]/p_ref

        grid_t_eq[:,:,k] .= (t_zero .- ΔT_y*sinθ2_2d  .- Δθ_z*cosθ2_2d.*log.(grid_p_norm)) .*  grid_p_norm.^kappa
        grid_t_eq[:,:,k] .= max.( t_strat, grid_t_eq[:,:,k])
    end
    
    for k=1:nd
        σ .= grid_p_full[:,:,k]./grid_ps

        #todo
        @assert(maximum(σ .- σ_b)/(1.0 - σ_b) <= 1.0)
        
        k_t .= k_a .+ (k_s - k_a)/(1.0-σ_b) * max.(0.0, σ .- σ_b) .* cosθ4_2d

        grid_δt_HS["HS_Newtonian"][:,:,k] .= -k_t .* (grid_t[:,:,k] - grid_t_eq[:,:,k])
        
    end
end 

!#######################################################################

function Rayleigh_Damping!(atmo_data::Atmo_Data, grid_p_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, 
    grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
    grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, 
    σ_b::Float64, k_f::Float64)
    
    #rayleigh damping of wind components near surface
    nλ, nθ, nd = size(grid_u)
    day_to_sec = 86400

    k_f = k_f/day_to_sec #s^-1

    grid_ps = grid_p_half[:,:,nd+1]
    
    
    σ = zeros(Float64, nλ, nθ)
    k_v = zeros(Float64, nλ, nθ)
    for k = 1:nd
        
        σ .= grid_p_full[:,:,k]./grid_ps

        #todo
        @assert(maximum(σ .- σ_b)/(1.0 - σ_b) <= 1.0)
        
        k_v = k_f/(1.0 - σ_b) * max.(0, σ .- σ_b)
        grid_δu[:,:,k]  .= -k_v .* grid_u[:,:,k]
        grid_δv[:,:,k]  .= -k_v .* grid_v[:,:,k]
        
    end
  
end 

!#######################################################################
"""Use linear response function (LRF) calc. from column rrtm"""
function Radiative_Forcing!(
    atmo_data::Atmo_Data, Δt::Int64, grid_p_full::Array{Float64,3}, 
    grid_t::Array{Float64,3}, grid_qv::Array{Float64,3}, 
    grid_δt_HS::Dict{String, Array{Float64, 3}};
    if_t_LRF::Bool=false, if_q_LRF::Bool=true
    )
    nλ, nθ, nd = size(grid_t)

    grid_δt_HS["LRF"] .*= 0
    file = load("/home/Quark/2025summer/prec/LRF/LRF_output_top300.dat")
    if if_q_LRF
        LRF_qv = file["LRF_SW_q"][:, :, :] .+ file["LRF_LW_q"][:, :, :] # run: 07/20 00:15
        # LRF_qv = file["LRF_LW_q"][:, :, :] # 7/21 19:19
        ref_qv = file["ref_q"][1, :, :]
        
        for ilat in 1:nθ
            LRF = LRF_qv[:, :, ilat]
            ref_qv_col = ref_qv[ilat, :]
            for ilon in 1:nλ
                pert_qv_col = grid_qv[ilon, ilat, :] - ref_qv_col
                grid_δt_HS["LRF"][ilon, ilat, :] .+= LRF * pert_qv_col
            end
        end
    end
    if if_t_LRF
        LRF_t = file["LRF_SW_T"][:, :, :] .+ file["LRF_LW_T"][:, :, :]
        ref_t = file["ref_T"][1, :, :]
        
        for ilat in 1:nθ
            LRF = LRF_t[:, :, ilat]
            ref_t_col = ref_t[ilat, :]
            for ilon in 1:nλ
                pert_t_col = grid_t[ilon, ilat, :] - ref_t_col
                grid_δt_HS["LRF"][ilon, ilat, :] .+= LRF * pert_t_col
            end
        end
    end

    grid_δt_HS["LRF"] ./= 86400 # k/day -> k/s
end
