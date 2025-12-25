using JGCM


function Atmos_Spectral_Dynamics_Main(physcis_params::Dict{String, Float64}, end_day::Int64 = 1200, spinup_day::Int64 = 200)
    # the decay of a sinusoidal disturbance to a zonally symmetric flow 
    # that resembles that found in the upper troposphere in Northern winter.
    # parameters
    name = "Spectral_Dynamics"
    num_fourier, nθ, nd = 21, 32, 10

    num_spherical = num_fourier + 1
    nλ = 2nθ
    
    radius = 6371000.0
    omega = 7.292e-5
    sea_level_ps_ref = 1.0e5
    init_t = 264.0

    ps_ref = sea_level_ps_ref
    t_ref = fill(300.0, nd)
    
    # Initialize mesh
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
    θc, λc = mesh.θc,  mesh.λc
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    
    vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
    # Initialize atmo_data
    do_mass_correction   = true
    do_energy_correction = true
    do_water_correction  = false
    
    use_virtual_temperature = false
    atmo_data = Atmo_Data(
        name, nλ, nθ, nd, 
        do_mass_correction, do_energy_correction, do_water_correction, 
        use_virtual_temperature, sinθ, radius, omega
    )
    
    # Initialize integrator
    damping_order = 4
    damping_coef = 1.15741e-4
    robert_coef  = 0.04 
    
    implicit_coef = 0.5
    day_to_sec = 86400
    start_time = 0
    end_time = end_day*day_to_sec  
    Δt = 1200
    ### CJY
    if warm_start_file_name != "None" 
        init_step = false # => In leapfrog would NOT do damping at initital time (should use in warm start case)
    else
        init_step = true # => In leapfrog would do damping at initital time
    end
    
    integrator = Filtered_Leapfrog(robert_coef, 
    damping_order, damping_coef, mesh.laplacian_eig,
    implicit_coef, Δt, init_step, start_time, end_time)

    semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data, integrator, ps_ref, t_ref, mesh.wave_numbers)
        
    # Initialize data
    dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)

    NT = Int64(end_time / Δt)
    Get_Topography!(dyn_data.grid_geopots)
    Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t, dyn_data.grid_geopots, dyn_data)
    
    # Initialize output manager
    op_man = Output_Manager(mesh, vert_coord, start_time, end_time, spinup_day)

    # first step
    Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data, physcis_params, L)
    Update_Init_Step!(semi_implicit)
    integrator.time += Δt
    Update_Output!(op_man, dyn_data, integrator.time)
    
    # main loop
    for i = 2:NT

        Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data, physcis_params)
        integrator.time += Δt
        Update_Output!(op_man, dyn_data, integrator.time)

    end

    return op_man
    
end

