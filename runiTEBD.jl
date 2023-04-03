using LinearAlgebra
using ITensors
using KrylovKit
using Dates
using QuadGK

include("canonical.jl")
include("miscellaneous.jl")
include("iMPS_functions.jl")
include("iTEBDmain.jl")


println("-----------------------------------------")
println(Dates.now())

β = 0.3
J = 1.0
h = 0.0

iTEBDmain(β, J, h)
