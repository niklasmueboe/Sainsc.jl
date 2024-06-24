module Utils

using Base.Broadcast: @__dot__
using PooledArrays: PooledArray
using Unzip: unzip

function categoricalcoordinates(x...)
    coordinates = PooledArray(collect(zip(x...)), Int32)
    return coordinates.refs, unzip(coordinates.pool)
end

stringcoordinates(x, y) = @. string(x) * "_" * string(y)
stringcoordinates(x...) = [join((string(j) for j in i), "_") for i in zip(x...)]

end # module Utils
