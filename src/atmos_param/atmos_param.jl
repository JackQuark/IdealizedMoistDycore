module atmos_param

    using JGCM: Atmo_Data, Semi_Implicit_Solver, Get_Δt

    include("HS_Forcing.jl")
    include("lscale_cond.jl")

    export HS_Forcing!
    export lscale_cond!

end