β = 0.1
J = 1.0

#=
                                    o   o   o   o
                                      X       X
                                    o   o   o   o
                                          X     
                                    o   o   o   o
                                      X       X
                                    o   o   o   o
=#
xmat = [exp(-β * (-J) * (s1*s2 + s2*s3 + s3*s4 + s4*s1)) for s1 in [-1,1], s2 in [-1,1], s3 in [-1,1], s4 in [-1,1]]

#=
             \      /
                 X
             /      \
- [λa] - [A] - [λb] - [B] - [λa] -
=#
d = 2
ai = Index(d, "a,site")
bi = Index(d, "b,site")

χ = 5
μ = Index(χ, "left, bond")
ν = Index(χ, "right, bond")

A = randomITensor(ai, μ, ν)
B = randomITensor(bi, μ, ν)
λA = randomITensor(μ, ν)
λB = randomITensor(μ, ν)

#=
    ul   ur
       X
    dl   dr
=#
ul = Index(d, "ul,site")
ur = Index(d, "ur, site")
dl = Index(d, "dl, site")
dr = Index(d, "dr, site")

X = ITensor(xmat, ul, ur, dl, dr)

function update1(A, B, λA, λB, X; kwargs...)
    ul = commonind(X, X, "ul,site")
    ur = commonind(X, X, "ur,site")
    dl = commonind(X, X, "dl,site")
    dr = commonind(X, X, "dr,site")
    ai = commonind(A, A, "a, site")
    bi = commonind(B, B, "b, site")

    it = ITensor(1.0)

    # λB A
    it *= prime(λβ,ν) * replaceind(A, μ, ν')

    # [λB A] λA
    prime!(it, ν)
    it *= replaceind(λA, μ, ν')

    # [λB A λA] B
    prime!(it, ν)
    it *= replaceind(B, μ, ν')

    # [λB A λA B] λB
    prime!(it, ν)
    it *= replaceind(λB, μ, ν')

    # [λB A λA B λB] X
    replaceinds!(it, [ai, bi], [dl, dr])
    it *= X

    #=
              ai bi
       μ [λB A λA B λB X] ν
    =#
    replaceinds!(it, [ul, ur], [ai, bi])

    # svd
    if :cutoff in keys(kwargs)
      cutoff = values(kwargs).cutoff
    else
      cutoff = 1e-8
    end
    
    U, S, V = svd(it, μ, ai; cutoff=1e-8)

    ν1 = commonind(U, S)
    replaceind!(U, ν1, ν)

    λBinv = ITensor(inv(array(λB, μ, ν)), μ, ν)

    λA = ITensor(S, μ, ν)

