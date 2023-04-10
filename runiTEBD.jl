using LinearAlgebra
using ITensors
using KrylovKit
using Dates
using QuadGK
using MAT

const βc = 0.5 * log(√2 + 1)

include("canonical.jl")
include("miscellaneous.jl")
include("iMPS_functions.jl")
include("iTEBDmain.jl")

println("-----------------------------------------")
println(Dates.now())
timeStamp = Dates.format(now(), "DyyyymmddTHHMMSS")

J = 1.0

lsβ = range(.94*βc, βc, length=64)
lsh = range(0, 0.01, length=64)

matfe = zeros(length(lsβ), length(lsh))

for (idxβ, β) in enumerate(lsβ)
    for (idxh, h) in enumerate(lsh)
        matfe[idxβ, idxh] = iTEBDmain(β, J, h; showQ = false)
    end
end

file = matopen(timeStamp*".mat", "w")
write(file, "matfe", matfe)
write(file, "lsb", collect(lsβ))
write(file, "lsh", collect(lsh))
close(file)
