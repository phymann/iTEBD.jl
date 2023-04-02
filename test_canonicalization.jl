#=
    Canonicalization of iMPS based on *PRB 78, 155117 (2008)*
=#

using ITensors
using KrylovKit
using LinearAlgebra
using PrettyTables
@ptconf tf = tf_borderless noheader = true crop = :horizontal formatters = ft_printf("%3.3e")

include("canonical.jl")
include("miscellaneous.jl")

function normalizeiMPS(Γ, λ)
   μ = commonind(Γ, Γ, tags="left, bond")
   ν = commonind(Γ, Γ, tags="right, bond")
   it1 = ITensor(1.)
   it1 *= replaceind(λ, ν, μ') * prime(Γ, "bond") * replaceind(λ, μ, ν')
   nrm = scalar(it1 * dag(it1))
   @show nrm
   Γ /= nrm^(1/2)
   # λ /= nrm^(1/6)

   μ = commonind(Γ, Γ, tags="left, bond")
   ν = commonind(Γ, Γ, tags="right, bond")
   it1 = ITensor(1.)
   it1 *= replaceind(λ, ν, μ') * prime(Γ, "bond") * replaceind(λ, μ, ν')
   nrm = scalar(it1 * dag(it1))
   if nrm ≈ 1
      println("now it is normalized!")
   else
      @infiltrate
   end
   return Γ, λ
end

let
   println("-----------------------------------------------------------------")
   println("random iMPS")
   d = 2
   χ = 12

   i = Index(d, "up,site")
   α = Index(χ, "left,bond")
   β = Index(χ, "right,bond")

   # ---------------------------------------------
   # initialize iMPS

   #=
            -α-[Γ]-β-
               |
               i
               |
   =#
   # Γ = randomITensor(ComplexF64, i, α, β)
   Γ = randomITensor(Float64, α, i, β)

   #       -α-[λ]-β-
   λ = diagITensor(rand(Float64, dim(α)), α, β)

   # println("--------------- first run ---------------")
   # check canonicalization and get R and L matrices
   R, L = canonQ(Γ, λ; showQ = true)
   Γ, λ = make_canon(Γ, λ, R, L)

   Γ, λ = normalizeiMPS(Γ, λ)

   # ---------------------------------------------
   # check canonicalization, again
   println("-----------------------------------------------------------------")
   println("After canonicalization process")
   R, L = canonQ(Γ, λ; showQ = true)

   return nothing
end