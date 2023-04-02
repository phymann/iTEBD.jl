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

#=
    ul   ur
       X
    dl   dr
=#
ul = Index(d, "ul,site")
ur = Index(d, "ur, site")
dl = Index(d, "dl, site")
dr = Index(d, "dr, site")

gateX = ITensor(xmat, ul, ur, dl, dr)

function update1(A, B, X)
    ul = commonind(X, X, "ul,site")
    ur = commonind(X, X, "ur,site")
    dl = commonind(X, X, "dl,site")
    dr = commonind(X, X, "dr,site")
