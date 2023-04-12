using LinearAlgebra
using ITensors
using KrylovKit
using Dates
using QuadGK
using MAT
using PrettyTables
@ptconf tf = tf_borderless noheader = true crop = :horizontal

const βc = 0.5 * log(√2 + 1)
const Tc = βc^-1

timeStamp = Dates.format(now(), "DyyyymmddTHHMMSS")
fname = "data/fe"*timeStamp*".mat"

dt = @elapsed let
    include("canonical.jl")
    include("miscellaneous.jl")
    include("iMPS_functions.jl")
    include("mainiTEBD.jl")

    println("-----------------------------------------")
    println(Dates.now())

    J = 1.0

    # nβ = 128
    # nh = 128

    # lsβ = range(.8*βc, 0.9995*βc, length=nβ)
    # lsh = range(0, 0.3, length=nh)

    ΔT = 0.005
    maxT = Tc+0.35
    minT = Tc+0.01
    lsT = minT:ΔT:maxT
    nT = length(lsT)
    lsβ = lsT.^-1
    nβ = nT
    @show nT

    Δh = 0.002
    maxh = 0.03
    minh = 0.001
    lsh = minh:Δh:maxh
    nh = length(lsh)
    @show nh

    matfe = zeros(nβ, nh)

    n_notconv = 0
    hls_notconv = []
    βls_notconv = []

    for (idxβ, β) in enumerate(lsβ)
        for (idxh, h) in enumerate(lsh)
            matfe[idxβ, idxh], convergenceQ = mainiTEBD(β, J, h; showQ=false, nrepeat=1024,
            maxdim=32, cutoff=1e-12)
            if !convergenceQ
                n_notconv += 1
                append!(hls_notconv, h)
                append!(βls_notconv, β)
            end
        end
    end

    file = matopen(fname, "w")
    write(file, "matfe", matfe)
    write(file, "lsb", collect(lsβ))
    write(file, "lsh", collect(lsh))
    write(file, "n_notconv", n_notconv)
    write(file, "hls_notconv", hls_notconv)
    write(file, "bls_notconv", βls_notconv)
    close(file)

    println("------------------------------")
    println("# of cases without convergence = $n_notconv")
    println("------------------------------")
    @pt hls_notconv
    println("------------------------------")
    @pt βls_notconv
end

file = matopen(fname, "r+")
write(file, "dt", dt)
close(file)

@show dt
@show fname
