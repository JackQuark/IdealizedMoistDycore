

function lscale_cond!(
    semi_implicit::Semi_Implicit_Solver,
    grid_p_full::Array{Float64, 3},
    grid_t_c::Array{Float64, 3}, 
    grid_q_c::Array{Float64, 3},
    grid_δt::Array{Float64, 3}, 
    grid_δq::Array{Float64, 3},
    factor3::Array{Float64, 3}, 
    L::Float64 = 0.2
)

    integrator = semi_implicit.integrator
    Δt         = Get_Δt(integrator)
    # These params are different to what we use in Atmo_Data. What a stupid...
    cp         = 1004. #
    Lv         = 2.5*10^6.
    Rd         = 287.04 
    Rv         = 461. #
    
    #grid_tracers_c_max     = deepcopy(grid_q_c)
    es                 = 611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t_c))
    grid_tracers_c_max = (0.622 .* es) ./ (grid_p_full .- 0.378 .* es) 
    dq_sat_dT          = Lv.*grid_tracers_c_max./ (Rv .*grid_t_c.^2)
    #@info "max: ", maximum(dq_sat_dT)
    # grid_d_full2 = dyn_data.grid_d_full2
    # @info maximum(grid_d_full2), minimum(grid_d_full2)

    ### Condensation_rate == grid_tracers_diff
    factor3  .= (max.(grid_q_c, grid_tracers_c_max) .- grid_tracers_c_max) ./ (1 .+ (Lv / cp) .* dq_sat_dT) ./(2 .* Δt)
    # grid_q_c       .-= (max.(grid_q_c, grid_tracers_c_max) .- grid_tracers_c_max) ./ (1 .+ (Lv / cp) .* dq_sat_dT)
    ### error ###
    grid_δq  .= copy(factor3)
    #############
    grid_δt  .= (factor3 .* Lv ./ cp) .* L 
    # latent heat feedback to temperature tendency 
end