# modified by Quark (summer2025)
export Output_Manager, Update_Output!, Finalize_Output!, All_to_Final

mutable struct Output_Manager
    nλ::Int64
    nθ::Int64
    nd::Int64
    n_day::Int64
    
    day_to_sec::Int64
    start_time::Int64
    end_time::Int64
    current_time::Int64
    spinup_day::Int64

    λc::Array{Float64, 1}
    θc::Array{Float64, 1}
    σc::Array{Float64, 1}

    op_data::Dict{String,Array}

    is1st_update::Bool
    final_keys::Vector{String} # required vars for warmstart file
end

function Output_Manager(
    mesh::Spectral_Spherical_Mesh, 
    vert_coord::Vert_Coordinate, 
    start_time::Int64, end_time::Int64, spinup_day::Int64, 
    num_grid_tracters::Int64=1, num_spe_tracters::Int64=1
) ### By CJY2
    
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    
    day_to_sec = 86400
    current_time = start_time

    λc = mesh.λc
    θc = mesh.θc

    #todo definition of sigma coordinate
    bk = vert_coord.bk
    σc = (bk[2:nd+1] + bk[1:nd])/2.0
  
    n_day = Int64((end_time - start_time)/ (day_to_sec/4) )

    grid_geopots_xyzt = zeros(Float64, nλ, nθ, 1, n_day)
    num_fourier, nθ, nd = 42, 64, 20
    num_spherical = num_fourier + 1
    #########################################################
    op_data = Dict{String,Array}(
        "spe_vor_c_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_vor_p_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_div_c_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_div_p_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_lnps_c_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day), ##
        "spe_lnps_p_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day), ##
        "spe_t_c_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_t_p_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_tracers_c_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "spe_tracers_p_xyzt" => zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day), ##
        "grid_u_c_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_u_p_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_v_c_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_v_p_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_ps_c_xyzt" => zeros(Float64, nλ, nθ, 1, n_day), ##
        "grid_ps_p_xyzt" => zeros(Float64, nλ, nθ, 1, n_day), ##
        "grid_t_c_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_t_p_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_tracers_n_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_tracers_c_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_tracers_p_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        # "grid_tracers_diff_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_δtracers_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "factor1_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "factor2_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "factor3_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "factor4_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_p_full_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), #
        "grid_p_half_xyzt" => zeros(Float64, nλ, nθ, nd+1, n_day), #
        "grid_geopots_xyzt" => zeros(Float64, nλ, nθ, 1, n_day), ##
        "grid_vor_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_div_xyzt" => zeros(Float64, nλ, nθ, nd, n_day), ##
        "grid_δu_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_δv_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_δt_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_z_full_xyzt" => zeros(Float64, nλ, nθ, nd, n_day),
        "grid_w_full_xyzt" => zeros(Float64, nλ, nθ, nd, n_day)
    )

    is1st_update = true

    final_keys = [
        "grid_geopots_xyzt", "grid_vor_xyzt", "grid_div_xyzt",
        "spe_vor_c_xyzt", "spe_vor_p_xyzt",
        "spe_div_c_xyzt", "spe_div_p_xyzt",
        "spe_lnps_c_xyzt", "spe_lnps_p_xyzt",
        "spe_t_c_xyzt", "spe_t_p_xyzt",
        "grid_u_c_xyzt", "grid_u_p_xyzt",
        "grid_v_c_xyzt", "grid_v_p_xyzt",
        "grid_ps_c_xyzt", "grid_ps_p_xyzt",
        "grid_t_c_xyzt", "grid_t_p_xyzt",
        "spe_tracers_c_xyzt", "spe_tracers_p_xyzt",
        "grid_tracers_n_xyzt", "grid_tracers_c_xyzt", "grid_tracers_p_xyzt"
    ]

    return Output_Manager(
        nλ, nθ, nd, n_day,
        day_to_sec, start_time, end_time, current_time, spinup_day, 
        λc, θc, σc,
        op_data, is1st_update,
        final_keys
    )
end

function Update_Output!(
    _op_man::Output_Manager, 
    dyn_data::Dyn_Data, 
    current_time::Int64
)
    @assert(current_time > _op_man.current_time)
    _op_man.current_time = current_time
    day_to_sec, start_time, n_day = _op_man.day_to_sec, _op_man.start_time, _op_man.n_day

    i_day = Int(div(current_time - start_time - 1, day_to_sec/4) + 1)

    if (i_day > n_day)
        @info "Warning: i_day > n_day in Output_Manager:Update!"
        return 
    end

    ############################################################
    # specral vor 
    _op_man.op_data["spe_vor_c_xyzt"][:,:,:,i_day] .= dyn_data.spe_vor_c[:,:,:]
    _op_man.op_data["spe_vor_p_xyzt"][:,:,:,i_day] .= dyn_data.spe_vor_p[:,:,:]
        
    # specral div
    _op_man.op_data["spe_div_c_xyzt"][:,:,:,i_day] .= dyn_data.spe_div_c[:,:,:]
    _op_man.op_data["spe_div_p_xyzt"][:,:,:,i_day] .= dyn_data.spe_div_p[:,:,:]
    
    # specral height or surface pressure
    _op_man.op_data["spe_lnps_c_xyzt"][:,:,:,i_day] .= dyn_data.spe_lnps_c[:,:,:]
    _op_man.op_data["spe_lnps_p_xyzt"][:,:,:,i_day] .= dyn_data.spe_lnps_p[:,:,:]

    # specral temperature
    _op_man.op_data["spe_t_c_xyzt"][:,:,:,i_day] .= dyn_data.spe_t_c[:,:,:]
    _op_man.op_data["spe_t_p_xyzt"][:,:,:,i_day] .= dyn_data.spe_t_p[:,:,:]
        
    # specral tracer
    _op_man.op_data["spe_tracers_c_xyzt"][:,:,:,i_day] .= dyn_data.spe_tracers_c[:,:,:]
    _op_man.op_data["spe_tracers_p_xyzt"][:,:,:,i_day] .= dyn_data.spe_tracers_p[:,:,:]
    ############################################################
    # grid w-e velocity
    # _op_man.op_data["grid_u_n_xyzt"][:,:,:,i_day]  .= dyn_data.grid_u_n[:,:,:] 
    _op_man.op_data["grid_u_c_xyzt"][:,:,:,i_day]  .= dyn_data.grid_u_c[:,:,:] 
    _op_man.op_data["grid_u_p_xyzt"][:,:,:,i_day]  .= dyn_data.grid_u_p[:,:,:]
    
    # grid n-s velocity
    # _op_man.op_data["grid_v_n_xyzt"][:,:,:,i_day]  .= dyn_data.grid_v_n[:,:,:]
    _op_man.op_data["grid_v_c_xyzt"][:,:,:,i_day]  .= dyn_data.grid_v_c[:,:,:]
    _op_man.op_data["grid_v_p_xyzt"][:,:,:,i_day]  .= dyn_data.grid_v_p[:,:,:]

    # grid surface pressure
    _op_man.op_data["grid_ps_c_xyzt"][:,:,:,i_day]  .= dyn_data.grid_ps_c[:,:,:]
    _op_man.op_data["grid_ps_p_xyzt"][:,:,:,i_day]  .= dyn_data.grid_ps_p[:,:,:]

    # grid temperature
    # _op_man.op_data["grid_t_n_xyzt"][:,:,:,i_day]  .= dyn_data.grid_t_n[:,:,:]
    _op_man.op_data["grid_t_c_xyzt"][:,:,:,i_day]  .= dyn_data.grid_t_c[:,:,:]
    _op_man.op_data["grid_t_p_xyzt"][:,:,:,i_day]  .= dyn_data.grid_t_p[:,:,:]

    # grid tracer
    _op_man.op_data["grid_tracers_n_xyzt"][:,:,:,i_day]  .= dyn_data.grid_tracers_n[:,:,:]
    _op_man.op_data["grid_tracers_c_xyzt"][:,:,:,i_day]  .= dyn_data.grid_tracers_c[:,:,:]
    _op_man.op_data["grid_tracers_p_xyzt"][:,:,:,i_day]  .= dyn_data.grid_tracers_p[:,:,:]

    # _op_man.op_data["grid_tracers_diff_xyzt"][:,:,:,i_day]  .= dyn_data.grid_tracers_diff[:,:,:]
    _op_man.op_data["grid_δtracers_xyzt"][:,:,:,i_day]  .= dyn_data.grid_δtracers[:,:,:]
    ############################################################
    # factors
    _op_man.op_data["factor1_xyzt"][:,:,:,i_day]  .= dyn_data.factor1[:,:,:]
    _op_man.op_data["factor2_xyzt"][:,:,:,i_day]  .= dyn_data.factor2[:,:,:]
    _op_man.op_data["factor3_xyzt"][:,:,:,i_day]  .= dyn_data.factor3[:,:,:]
    _op_man.op_data["factor4_xyzt"][:,:,:,i_day]  .= dyn_data.factor4[:,:,:]
    
    # pressure
    _op_man.op_data["grid_p_full_xyzt"][:,:,:,i_day]  .= dyn_data.grid_p_full[:,:,:]
    _op_man.op_data["grid_p_half_xyzt"][:,:,:,i_day]  .= dyn_data.grid_p_half[:,:,:]
    
    # geopotential
    _op_man.op_data["grid_geopots_xyzt"][:,:,:,i_day]  .= dyn_data.grid_geopots[:,:,:]
    ############################################################
    
    _op_man.op_data["grid_vor_xyzt"][:,:,:,i_day]  .= dyn_data.grid_vor[:,:,:]
    _op_man.op_data["grid_div_xyzt"][:,:,:,i_day]  .= dyn_data.grid_div[:,:,:]

    # Tendency
    _op_man.op_data["grid_δu_xyzt"][:,:,:,i_day]  .= dyn_data.grid_δu[:,:,:]
    _op_man.op_data["grid_δv_xyzt"][:,:,:,i_day]  .= dyn_data.grid_δv[:,:,:]
    _op_man.op_data["grid_δt_xyzt"][:,:,:,i_day]  .= dyn_data.grid_δt[:,:,:]

    _op_man.op_data["grid_z_full_xyzt"][:,:,:,i_day]  .= dyn_data.grid_z_full[:,:,:]
    _op_man.op_data["grid_w_full_xyzt"][:,:,:,i_day]  .= dyn_data.grid_w_full[:,:,:]

    # HS forcing var keys generate during first update
    if _op_man.is1st_update
        for k in keys(dyn_data.grid_δt_HS)
            _op_man.op_data[k] = zeros(Float64, _op_man.nλ, _op_man.nθ, _op_man.nd, n_day)    
        end
        is1st_update = false
    end
    # load forcing from SH_forcing.jl (see SH_forcing.jl for details)
    for k in keys(dyn_data.grid_δt_HS)
        _op_man.op_data[k][:,:,:,i_day] .= dyn_data.grid_δt_HS[k][:,:,:]
    end

end

function Finalize_Output!(
    _op_man::Output_Manager, 
    ofname_final::String = "None", ofname_all::String = "None"
)
    n_day = _op_man.n_day
    op_data = _op_man.op_data
    final_keys = _op_man.final_keys

    if ofname_all != "None"
        JLD2.save(ofname_all, op_data)  
    end

    if ofname_final != "None"
        final_data = Dict{String, Any}()
        for k in final_keys
            if haskey(op_data, k)
                final_k = replace(k, "_xyzt" => "_final")
                final_data[final_k] = @view op_data[k][:, :, :, end]
            end
        end
        JLD2.save(ofname_final, final_data)
    end
end

"""
    All_to_Final(all_fname::String, final_fname::String)
    
Extract the final time snapshot from an "all times" JLD2 file and save it to a new JLD2 file.
```
"""
function All_to_Final(
    all_fpath::String, final_fpath::String
) # add by Quark
    f = load(all_fpath)
    fin_data = Dict{String,Array}(
        "grid_geopots_final" => f["grid_geopots_xyzt"][:,:,:,end],
        "grid_vor_final"  => f["grid_vor_xyzt"][:,:,:,end],
        "grid_div_final"  => f["grid_div_xyzt"][:,:,:,end],
        "spe_vor_c_final"    => f["spe_vor_c_xyzt"][:,:,:,end],
        "spe_vor_p_final"    => f["spe_vor_p_xyzt"][:,:,:,end],
        "spe_div_c_final"    => f["spe_div_c_xyzt"][:,:,:,end],
        "spe_div_p_final"    => f["spe_div_p_xyzt"][:,:,:,end],
        "spe_lnps_c_final"   => f["spe_lnps_c_xyzt"][:,:,:,end],
        "spe_lnps_p_final"   => f["spe_lnps_p_xyzt"][:,:,:,end],
        "spe_t_c_final"      => f["spe_t_c_xyzt"][:,:,:,end],
        "spe_t_p_final"      => f["spe_t_p_xyzt"][:,:,:,end],
        "grid_u_c_final"  => f["grid_u_c_xyzt"][:,:,:,end],
        "grid_u_p_final"  => f["grid_u_p_xyzt"][:,:,:,end],
        "grid_v_c_final"  => f["grid_v_c_xyzt"][:,:,:,end],
        "grid_v_p_final"  => f["grid_v_p_xyzt"][:,:,:,end],
        "grid_ps_c_final" => f["grid_ps_c_xyzt"][:,:,:,end],
        "grid_ps_p_final" => f["grid_ps_p_xyzt"][:,:,:,end],
        "grid_t_c_final"  => f["grid_t_c_xyzt"][:,:,:,end],
        "grid_t_p_final"  => f["grid_t_p_xyzt"][:,:,:,end],
        "spe_tracers_c_final"  => f["spe_tracers_c_xyzt"][:,:,:,end],
        "spe_tracers_p_final"  => f["spe_tracers_p_xyzt"][:,:,:,end],
        "grid_tracers_n_final" => f["grid_tracers_n_xyzt"][:,:,:,end],
        "grid_tracers_c_final" => f["grid_tracers_c_xyzt"][:,:,:,end],
        "grid_tracers_p_final" => f["grid_tracers_p_xyzt"][:,:,:,end],
    )
    f.close()
    JLD2.save(final_fpath, fin_data)
    println("Saving final snapshot to: ", final_fpath)
end
