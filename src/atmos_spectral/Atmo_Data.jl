export Atmo_Data, Compute_Abs_Vor!


struct Atmo_Data
    name::String                    # unique identifier

    nλ::Int64                       # zonal grids
    nθ::Int64                       # meridional grids
    nd::Int64                       # vertical grids
    
    do_mass_correction::Bool
    do_energy_correction::Bool
    do_water_correction::Bool
    use_virtual_temperature::Bool

    radius::Float64                 # planet radius
    omega::Float64                  # planet angular speed
    grav::Float64                   # gravity
    rdgas::Float64                  # dry ideal gas constant
    kappa::Float64                  # kappa
    rvgas::Float64                  # moist ideal gas constant
    cp_air::Float64                 # Isobaric specific heat of air
    Lv::Float64                     # Latent heat of vaporization

    coriolis::Array{Float64,1}      # Coriolis

end

function Atmo_Data(name::String,  
                   nλ::Int64, nθ::Int64, nd::Int64, 
                   do_mass_correction::Bool, do_energy_correction::Bool, do_water_correction::Bool, use_virtual_temperature::Bool,
                   sinθ::Array{Float64,1},
                   radius::Float64=6.371e6, 
                   omega::Float64=7.292e-5, grav::Float64=9.80, kappa::Float64=2.0/7.0, 
                   rdgas::Float64=287.04, rvgas::Float64=461.50)
    """
    Constructs a new `Atmo_Data` object,
    initializing all atmospheric configuration parameters, logic flags, physical constants, and derived dynamical fields.
    """
    coriolis = 2 * omega * sinθ
    cp_air = rdgas/kappa
    Lv = 2.5e6
    
    Atmo_Data(name,  
              nλ, nθ, nd,
              do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature,
              radius, omega, grav, rdgas, kappa, rvgas, cp_air, Lv,
              coriolis)
end



function Compute_Abs_Vor!(grid_vor::Array{Float64,3}, coriolis::Array{Float64,1}, grid_absvor::Array{Float64,3})
    nλ, nθ, nd = size(grid_vor)

    for j=1:nθ
        grid_absvor[:,j,:] .= grid_vor[:,j,:] .+ coriolis[j]
    end
end