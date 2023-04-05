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

J = 1.0

lsβ = 0.2:0.005:βc
lsh = -0.5:0.005:0.5

matfe = zeros(length(lsβ), length(lsh))

for (idxβ, β) in enumerate(lsβ)
    for (idxh, h) in enumerate(lsh)
        matfe[idxβ, idxh] = iTEBDmain(β, J, h; showQ = false)
    end
end

file = matopen("matfile.mat", "w")
write(file, "matfe", matfe)
write(file, "lsb", collect(lsβ))
write(file, "lsh", collect(lsh))
close(file)
