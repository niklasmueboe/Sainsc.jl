module SingleCellProjectionsExt

using SingleCellProjections: DataMatrix

using Sainsc

function Sainsc.getlocalmaxima(::Type{DataMatrix}, counts, localmax, kernel; genes=nothing)
    return DataMatrix(getlocalmaxima(counts, localmax, kernel; genes=genes))
end

function Sainsc.readstereoseqbinned(::Type{DataMatrix}, file, binsize::Integer)
    return DataMatrix(readstereoseqbinned(file, binsize))
end

end # module SingleCellProjectionsExt
