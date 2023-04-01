using LinearAlgebra
using ITensors
using PrettyTables
using KrylovKit
using Infiltrator

include("canonical.jl")
include("miscellaneous.jl")


"""
To update Γ and λ
"""
function updateit(Γ, λ, a)

  μ = commonind(Γ, Γ, tags="left, bond")
  ν = commonind(Γ, Γ, tags="right, bond")
  j = commonind(a, a, tags="left, site")
  l = commonind(a, a, tags="right, site")

  it = ITensor(1.)
  #=
  --------------------------------------------------
  μ-[λΓλ]-ν
    |
    k
  =#
  it *= replaceind(λ, ν, μ') * prime(Γ, "bond") * replaceind(λ, μ, ν')

  #=
  --------------------------------------------------
  μ-[λΓλ]-ν
    |
  j-[a]-l
    |
    k
  =#
  it *= a

  #=
  --------------------------------------------------
  μ-[λΓλ]-ν
  [ a ]-l
    |
    k
  =#
  CL = combiner(j, μ)
  cli = combinedind(CL)
  it *= CL
  settags!(it, "left, bond", cli)

  λ *= delta(j, l)
  λ *= CL
  settags!(λ, "left, bond", cli)

  #=
  --------------------------------------------------
  μ-[λΓλ]-ν
  [ a ]
    |
    k
  =#
  CR = combiner(l, ν)
  cri = combinedind(CR)
  it *= CR
  settags!(it, "right, bond", cri)

  λ *= CR
  settags!(λ, "right, bond", cri)

  #=
  --------------------------------------------------
  μ-[λΓλ]-ν
  [ a ]
    |
    i
  =#
  replacetags(it, "down, site", "up, site")

  # Γ, λ = normalizeiMPS(Γ, λ, μ, ν)

  return it, λ
end

"""

"""
function getZ(Γ, λ, a)
  μ = commonind(Γ, Γ, tags="left, bond")
  ν = commonind(Γ, Γ, tags="right, bond")
  l = commonind(a, a, tags="right, site")
  j = commonind(a, a, tags="left, site")
  it1 = ITensor(1.)
  sqrtλ = ITensor(sqrt(array(λ, μ, ν)), μ, ν)
  it1 *= replaceind(sqrtλ, ν, μ') * prime(Γ, "bond") * replaceind(sqrtλ, μ, ν')
  it2 = dag(prime(it1, "bond"))
  replaceind!(it2, i, k)
  vals, _, _, = eigsolve(x -> array(((ITensor(x, ν, l, ν') * it1) * a) * it2, μ, j, μ'), rand(dim(ν), dim(l), dim(ν')))
  return vals[1]
end


"""
find norm of iMPS
"""
function normalizeiMPS(Γ, λ)
  μ = commonind(Γ, Γ, tags="left, bond")
  ν = commonind(Γ, Γ, tags="right, bond")
  it1 = ITensor(1.)
  it1 *= prime(λ, ν) * prime(Γ, "bond") * prime(λ, μ)
  nrm = scalar(it1 * dag(it1))
  Γ /= nrm^(1/6)
  λ /= nrm^(1/6)
  return Γ, λ
end




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

    In order to obtain Z,
    we will first find {Γ, λ} for Ψu and Ψd, by iteratively applying
    the transfer matrix T on an initial state Ψ0
    from Γ, λ, we can construct W
    then obtain the dominant eigenvalue ω of W
=#

let
  β = 0.1
  J = 1

  d = 2
  i = Index(d, "up, site")
  j = Index(d, "left, site")
  k = Index(d, "down, site")
  l = Index(d, "right, site")

  χ = 5
  μ = Index(χ, "left, bond")
  ν = Index(χ, "right, bond")

  #=
      μ-[Γ]-ν
          |
          i
  =#
  Γ = randomITensor(Float64, μ, i, ν)
  # ones()
  # Γ = ITensor()
  #=
        μ-[λ]-ν
  =#

  # λ = diagITensor(diag(array(randomITensor(Float64, μ, ν))), μ, ν)
  λ = ITensor(Diagonal(ones(χ)), μ, ν)
  R, L = canonQ(Γ, λ; showQ = true)
  Γ, λ = make_canon(Γ, λ, R, L)
  Γ, λ = normalizeiMPS(Γ, λ)
  canonQ(Γ, λ; showQ = true)

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

  for counti in 1:100
    @show counti
    Γold = Γ
    λold = λ
    Γ, λ = updateit(Γold, λold, a)
    R, L, torf = canonQ(Γ, λ; showQ = true)
    if ~torf
      println("update to a non-canonical iMPS")
    end
    Γ, λ = make_canon(Γ, λ, R, L; maxdim = χ, cutoff = 1e-2)


    # Γ, λ = normalizeiMPS(Γ, λ, μ, ν)
    @show norm(λ)
    @show norm(λold)
    cond = norm(λ - λold)/norm(λ)
    @show cond
    if cond < 1e-3
      break
    end
  end

  Γ, λ = normalizeiMPS(Γ, λ)

  return getZ(Γ, λ, a)
end