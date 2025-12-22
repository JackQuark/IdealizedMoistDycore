using JGCM

include("HS.jl")

end_day = 1
spinup_day = 0

L = 0.2
warm_start_file_name = "/data92/Quark/ctrl_2000d/startfrom_475day_final.dat"
initial_day = 0 # warm start day

physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 60.0, "Δθ_z" => 10.0) 
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day, L)

Finalize_Output!(op_man, "test_final.h5", "test_all.h5")

