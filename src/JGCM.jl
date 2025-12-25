__precompile__(true)
module JGCM

using Revise
using LinearAlgebra
using FFTW
using Statistics 
using JLD2
using MAT
using Interpolations
using TOML

include("atmos_spectral/Atmo_Data.jl")
include("atmos_spectral/Dyn_Data.jl")
include("atmos_spectral/Gauss_And_Legendre.jl")
include("atmos_spectral/Spectral_Spherical_Mesh.jl")
include("atmos_spectral/Vert_Coordinate.jl")
include("atmos_spectral/Press_And_Geopot.jl")
include("atmos_spectral/Time_Integrator.jl")
include("atmos_spectral/Semi_Implicit.jl") 

# shared
include("shared/Output_Manager.jl")

# Params
include("atmos_param/atmos_param.jl")
using .atmos_param

# Applications
include("atmos_spectral_barotropic/Barotropic_Dynamics.jl")
include("atmos_spectral_shallow/Shallow_Water_Dynamics.jl")
include("atmos_spectral/Spectral_Dynamics.jl")
end