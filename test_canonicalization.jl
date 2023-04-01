#=
    Canonicalization of iMPS based on *PRB 78, 155117 (2008)*
=#

using ITensors
using KrylovKit
using LinearAlgebra
using PrettyTables
@ptconf tf = tf_borderless noheader = true crop = :horizontal formatters = ft_printf("%3.3e")

function jwstop()
   println("---------------------------")
   println("---------------------------")
   println("---------------------------")
   println("--------No Problem!--------")
   println("---------------------------")
   println("---------------------------")
   println("---------------------------")
   error("stop here")
end

function jwdisplay(x,str)
   println("---------------------------")
   println("$str is ")
   display(x)
end

function jwchk(x::Bool)
   if !x
       error("sanity check failed")
   end
end

function canonQ(Γ, λ, α, β; showQ = false)
  #=
   contruct R matrix
               ------
                  |
                 [R]
                  |
               ------

                  ||

               -α-[Γ-λ]-β-
                  |
            -α'-[Γ*-λ*]-β-
   =#
   ε = 1E-8
   R_tmp = prime(Γ, β) * replaceind(λ, α, β')         # i, α, β
   R_tmp1 = dag(replaceinds(R_tmp, [α, β], [α', β'])) # i, α', β'
   R = R_tmp * R_tmp1                                 # α, α', β, β'

   # check whether iMPS is right canonical
   R1 = R * delta(β, β')
   R1mat = array(R1, α, α')
   R1mat0 = Diagonal([R1mat[1,1] for _ in 1:size(R1mat)[1]])
   if norm(R1mat - R1mat0) > ε
      println("--------------------------------")
      println("η*Id is NOT R's right eigenvector!")
      if showQ
         @pt R1mat
      end
   else
      println("⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ")
      println("η*Id is R's right eigenvector!")
      if showQ
         @pt R1mat
      end
   end

   #=
   construct L matrix
                     ------
                        |
                       [L]
                        |
                     ------

                        ||

                     -α-[λ-Γ]-β-
                         |
                         i
                         |
                    -α'-[λ-Γ]-β'-
   =#
   L_tmp = prime(Γ, α) * replaceind(λ, β, α')         # i, α, β
   L_tmp1 = dag(replaceinds(L_tmp, [α, β], [α', β'])) # i, α', β'
   L = L_tmp * L_tmp1                                 # α, α', β, β'

   #check whether iMPS is canonical:
   L1 = L * delta(α, α')
   L1mat = array(L1, β, β')
   L1mat0 = Diagonal([L1mat[1,1] for _ in 1:size(L1mat)[1]])
   if norm(L1mat - L1mat0) > ε
      println("-------------------------------")
      println("η*Id is NOT L's left eigenvector!")
      if showQ
         @pt L1mat
      end
   else
      println("⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ")
      println("η*Id is L's left eigenvector!")
      if showQ
         @pt L1mat
      end
   end

   return R, L
end

function make_canon(Γ, λ, R, L, α, β)
   ## ====================
   ## step 1: find X and Y
   ## ====================

   # ---------------------
   ### find X
   #### find the dominant right eigenvector (using anonymous function as input)
   vals, vecs, _ = eigsolve(x -> array(ITensor(x, β, β') * R, α, α'), rand(dim(β), dim(β')))
   η = vals[1]
   Vᵣvec = vecs[1]
   ##### making the vector real, if possible
   idx = findfirst(==(maximum(abs.(Vᵣvec))), abs.(Vᵣvec))
   Vᵣvec *= conj(Vᵣvec[idx])/abs(Vᵣvec[idx])

   #### switch back to ITensor
   Vᵣ = ITensor(real(Vᵣvec), β, β') # Vᵣ: β, β'
   jwchk(norm(R * Vᵣ - η * replaceinds(Vᵣ, [β, β'] ,[α, α'])) < 1e-8 * norm(Vᵣ))

   #### then do the decomposition
   # Vᵣ * W = W * D
   D, W = eigen(Vᵣ, β, β', ishermitian=true)
   comidx = commonind(W, D)   # comidx
                              # W: β', comidx
   uniqidxX = uniqueind(D, W) # D: comidx, uniqidxX
   Wl = replaceinds(W, [β', comidx], [β, uniqidxX])
   Dmat = array(D, comidx, uniqidxX)
   sqrtD = ITensor(sqrt(Dmat), comidx, uniqidxX) # sqrtD: comidx, uniqidxX
   # X = W * √D
   X = replaceind(W, β', β) * sqrtD # X:  β, uniqidxX
   Xp = dag(replaceind(X, β, β'))   # Xp: β', uniqidxX
   # X * X† = Vᵣ
   jwchk(norm(X*Xp - Vᵣ) < 1e-8 * norm(Vᵣ))

   # ---------------------
   ### find Y
   #### find the dominant *left* eigenvector (using anonymous function as input)
   vals, vecs, _ = eigsolve(x -> array(ITensor(x, α, α') * L, β, β'), rand(dim(α), dim(α')))
   jwchk(vals[1] ≈ η)
   Vₗvec = vecs[1]
   ##### making the vector real, if possible
   idx = findfirst(==(maximum(abs.(Vₗvec))), abs.(Vₗvec))
   Vₗvec *= conj(Vₗvec[idx])/abs(Vₗvec[idx])

   #### switch back to ITensor
   Vₗ = ITensor(Vₗvec, α, α')  # Vₗ: α, α'
   jwchk(norm(L * Vₗ - η * replaceinds(Vₗ, [α, α'], [β, β'])) < 1e-8 * norm(Vₗ))

   #### then do the decomposition
   # Vₗ * W = W * D
   D, W = eigen(Vₗ, α, α')
   comidx = commonind(W, D)   # comidx
                              # W: α', comidx
   uniqidxY = uniqueind(D, W) # D: comidx, uniqidxY
   Dmat = array(D, comidx, uniqidxY)
   sqrtD = ITensor(sqrt(Dmat), comidx, uniqidxY)
   # Y = W * √D
   Y = replaceind(W, α', α) * sqrtD # Y:  α, uniqidxY
   Yp = dag(replaceind(Y, α, α'))   # Yp: α', uniqidxY
   # Y * Y† = Vₗ
   jwchk(Y * Yp ≈ Vₗ)

   ## ====================================
   ## step 2: SVD => transpose(Y)λX = Uλ'V
   ## ====================================
   YtλX = replaceind(Y, α', α) * λ * X # YtλX: uniqidxY, uniqidxX
   U, S, V = svd(YtλX, uniqidxY)
   jwchk(norm(U * S * V - YtλX) < 1e-8 * norm(YtλX))
   αp = commonind(U, S) # U: uniqidxY, αp
   βp = commonind(V, S) # V: βp, uniqidxX
   λ = replaceinds(S, [αp, βp], [α, β])

   ## ======================================
   ## step 3: construct Γ': V invX Γ invYt U
   ## ======================================
   ### be very careful with the ordering of indices!!!
   Xinv = ITensor(inv(array(X, β, uniqidxX)), uniqidxX, α)
   Ytinv = ITensor(inv(array(Y, uniqidxY, α)), β, uniqidxY)
   Γ = V * Xinv * Γ * Ytinv * U
   ### be very careful with the ordering of indices!!!
   replaceinds!(Γ, [αp, βp], [β, α])

   return Γ, λ
end

function normalizeiMPS(Γ, λ, μ, ν)
   it1 = ITensor(1.)
   it1 *= prime(λ, ν) * prime(Γ, "bond") * prime(λ, μ)
   nrm = scalar(it1 * dag(it1))
   Γ /= (nrm^(1/6))
   λ /= nrm^(1/6)
   return Γ, λ
 end

let
   println("-----------------------------------------------------------------")
   println("random iMPS")
   d = 2
   χ = 12

   i = Index(d, "i,site")
   α = Index(χ, "α,bond")
   β = Index(χ, "β,bond")

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
   R, L = canonQ(Γ, λ, α, β; showQ = true)
   Γ, λ = make_canon(Γ, λ, R, L, α, β)

   Γ, λ = normalizeiMPS(Γ, λ, α, β)

   # ---------------------------------------------
   # check canonicalization, again
   println("-----------------------------------------------------------------")
   println("After canonicalization process")
   R, L = canonQ(Γ, λ, α, β; showQ = true)



   return nothing
end