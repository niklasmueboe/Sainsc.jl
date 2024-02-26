module SingleCellProjectionsExt

using SingleCellProjections: DataMatrix

using StereoSSAM

function StereoSSAM.getlocalmaxima(
    T::Type{DataMatrix}, counts, localmax, kernel; genes=nothing
)
    return DataMatrix(getlocalmaxima(counts, localmax, kernel; genes=genes))
end

function StereoSSAM.readstereoseqbinned(T::Type{DataMatrix}, file, s::Integer)
    return DataMatrix(readstereoseqbinned(file, s))
end

end # module SingleCellProjectionsExt
