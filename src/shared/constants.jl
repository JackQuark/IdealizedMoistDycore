module AtmosConstants

    export cnst # public constant object

    abstract type AbstractConstants end

    # default value sets =====( 
        # @kwdef for default values & flexible Float Type (32,64)
        Base.@kwdef struct CJYConstants{T<:AbstractFloat} <: AbstractConstants
            # refer to https://github.com/PeterChang1023/Moist_Dycore/tree/For_other_user
            radius::T = 6371.0e3
            omega::T = 7.292e-5
            g::T = 9.8
            p0::T = 100000.0
            rdgas::T = 287.04
            rvgas::T = 461.50
            kappa::T = 2.0 / 7.0
            cp::T = rdgas / kappa
        end

        Base.@kwdef struct IscaConstants{T<:AbstractFloat} <: AbstractConstants
            # TODO: https://github.com/ExeClim/Isca/blob/master/src/shared/constants/constants.F90
            ...
        end
    # )===== default value sets 

    # Load toml setting
    # TODO: config should be loaded at model.jl
    config = TOML.parsefile(WORKDIR + "/config.toml")
    
    if haskey(config, "constants")
        kwargs = Dict{Symbol, float_type}()
        
        for (k, v) in config["constants"]
            kwargs[Symbol(k)] = convert(float_type, v) # convert to mapping[key::Symbol -> value::float_type]
        end
        
        # TODO: allow quick selection of default value sets by toml
        const cnst = CJYConstants{Float64}(; kwargs...)
    else
        const cnst = CJYConstants{Float64}()
    end
    
end