using LinearAlgebra
using ITensors
using PrettyTables

include("canonical.jl")


let
#=
    two-dimensional Ising model
        H = ∑ h

        Z(β) = ∑ ∏ exp(-β h)

    define the Hermitian matrix Q
        Q[s,s'] = exp(-β h)
    for each bond of the original lattice

      |   |   |
    - o - o - o -
      |   |   |
    - o - o - o -
      |   |   |
    - o - o - o -
      |   |   |

    for each site, we define a order-4 tensor a
        a = ∑ₛ √Q(i,s) √Q(j,s) √Q(k,s) √Q(l,s)

    The partition function = contraction of an infinite 2D tensor network specified by a single tensor a

    Z(β) = lim(p→∞) tr(Tᵖ) = lim(p→∞) θᵖ
    θ = tr(W^q) = lim(q→∞) ω^q

    Z = lim(p,q→∞) ω^{pq}
=#

β = 0.1

d = 2
i = Index(d, "i, site")
j = Index(d, "j, site")
k = Index(d, "k, site")
l = Index(d, "l, site")

χ = 20
μ = Index(χ, "μ, bond")
ν = Index(χ, "ν, bond")

#=
     μ-[Γ]-ν
        |
        i
=#
Γ = randomITensor(Float64, μ, i, ν)
#=
      μ-[λ]-ν
=#
λ = diagITensor(rand(Float64, dim(μ)), μ, ν)

R, L = canonQ(Γ, λ, α, β; showQ = true)
Γ, λ = make_canon(Γ, λ, R, L, μ, ν)

Qmat = [ exp(-β * (-J)) exp(β * J); 
         exp(β * J) exp(-β * (-J))
    ]

Qi = ITensor(Qmat, i, i')
Qj = ITensor(Qmat, j, j')
Qk = ITensor(Qmat, k, k')
Ql = ITensor(Qmat, l, l')

delta4 = delta(i', j', k', l')

#=
      i
      |
   j-[a]-l
      |
      k
=#
a = Qi * Qj * Qk * Ql * delta4

#=
    To update 
=#
function updateit(Γ, λ, a)
  it = ITensor(1.)
  #=
  μ-[λΓλ]-ν
      |
      i
  =#
  it *= prime(λ, β) * prime(Γ, "bond") * prime(λ, α)
  #=
  μ-[λΓλ]-ν
      |
   j-[a]-l
      |
      k
  =#
  it *= a
  #=
  μ-[λΓλ]-ν
    [ a ]-l
      |
      k
  =#
  CL = combiner(j, μ)
  cli = combinedind(CL)
  it *= CL
  settags(it, "μ, bond"; tags=tags(cli))
  #=
  μ-[λΓλ]-ν
    [ a ]
      |
      k
  =#
  CR = combiner(l, ν)
  cri = combinedind(CR)
  it *= CR
  settags(it, "ν, bond"; tags=tags(cri))
  #=
  μ-[λΓλ]-ν
    [ a ]
      |
      i
  =#
  replacetags(it, "k, site", "i, site")
  return Γ, λ
end

function getZ(Γ, λ, Γold, λold, a)
  
end

for _ in 1:100
  Γold = Γ
  λold = λ
  Γ, λ = updateit(Γold, λold, a)
  R, L = canonQ(Γ, λ, α, β; showQ = false)
  Γ, λ = make_canon(Γ, λ, R, L, μ, ν; maxdim = dim(χ))


end