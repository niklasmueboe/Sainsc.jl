module AnnDataExt

using Muon: AnnData

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{AnnData}, counts, localmax, kernel; genes=nothing
)
    X, genes, obs = getlocalmaxima(counts, localmax, kernel; genes=genes)

    return AnnData(;
        X=permutedims(X),
        var_names=genes.gene,
        obs_names=obs.id,
        obsm=Dict("spatial" => hcat(obs.x, obs.y)),
    )
end

function StereoSSAM.readstereoseqbinned(T::Type{AnnData}, file, s::Integer)
    X, genes, obs = readstereoseqbinned(file, s)

    return AnnData(;
        X=permutedims(X),
        var_names=genes.gene,
        obs_names=obs.bin_id,
        obsm=Dict("spatial" => hcat(obs.x, obs.y)),
    )
end

end # module AnnDataExt
