using LinearAlgebra
using ITensors
using KrylovKit
using Dates
using QuadGK
using MAT
using PrettyTables
@ptconf tf = tf_borderless noheader = true crop = :horizontal

const βc = 0.5 * log(√2 + 1)

dt = @elapsed let
include("canonical.jl")
include("miscellaneous.jl")
include("iMPS_functions.jl")
include("mainiTEBD.jl")

    println("-----------------------------------------")
    println(Dates.now())
    timeStamp = Dates.format(now(), "DyyyymmddTHHMMSS")

    J = 1.0

    nβ = 16
    nh = 128

    lsβ = range(.9*βc, βc, length=nβ)
    lsh = range(0, 0.01, length=nh)

    matfe = zeros(nβ, nh)

    n_notconv = 0
    hls_notconv = []
    βls_notconv = []

    for (idxβ, β) in enumerate(lsβ)
        for (idxh, h) in enumerate(lsh)
            matfe[idxβ, idxh], convergenceQ = mainiTEBD(β, J, h; showQ=false, nrepeat=1024)
            if !convergenceQ
                n_notconv += 1
                append!(hls_notconv, h)
                append!(βls_notconv, β)
            end
        end
    end

    println("------------------------------")
    println("# of cases without convergence = $n_notconv")
    println("------------------------------")
    @pt hls_notconv
    println("------------------------------")
    @pt βls_notconv
end

file = matopen(timeStamp*".mat", "w")
write(file, "matfe", matfe)
write(file, "lsb", collect(lsβ))
write(file, "lsh", collect(lsh))
write(file, "n_notconv", n_notconv)
write(file, "hls_notconv", hls_notconv)
write(file, "bls_notconv", βls_notconv)
write(file, "dt", dt)
close(file)
