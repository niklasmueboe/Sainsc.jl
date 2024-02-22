module AnnDataExt


export getlocalmaximaanndata


using Muon: AnnData

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{AnnData},
    sp_c,
    localmax,
    kernel;
    genes = nothing,
)
    X, genes, obs_names, coordinates =
        StereoSSAM._getlocalmaxima(sp_c, localmax, kernel; genes = genes)

    AnnData(
        X = X,
        var_names = genes,
        obs_names = obs_names,
        obsm = Dict("spatial" => hcat(coordinates...)),
    )
end

end # module AnnDataExt
