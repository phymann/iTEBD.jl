function canonQ(Γ, λ; showQ = false, checkQ = false)
    # parameters
    ε = 1E-3

    # setup indices
    α = commonind(Γ, Γ, tags="left, bond")
    β = commonind(Γ, Γ, tags="right, bond")

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
        -α'-[Γ*-λ]-β-
    =#
    R_tmp = prime(Γ, β) * replaceind(λ, α, β')         # i, α, β
    R_tmp1 = dag(replaceinds(R_tmp, [α, β], [α', β'])) # i, α', β'
    R = R_tmp * R_tmp1                                 # α, α', β, β'

    # check whether iMPS is right canonical
    R1 = R * delta(β, β')
    R1mat = array(R1, α, α')
    R1mat0 = Diagonal([R1mat[1,1] for _ in 1:size(R1mat)[1]])

    torf = true
    if norm(R1mat - R1mat0)/norm(R1mat) > ε
        torf = false
        # @show norm(R1mat - R1mat0)
        if showQ
            println("--------------------------------")
            println("Id is NOT R's right eigenvector!")
            @show norm(R1mat - R1mat0)/norm(R1mat)
        end
    else
        if showQ
            println("⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ")
            println("Id is R's right eigenvector!")
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
                    -α'-[λ-Γ*]-β'-
    =#
    L_tmp = prime(Γ, α) * replaceind(λ, β, α')         # i, α, β
    L_tmp1 = dag(replaceinds(L_tmp, [α, β], [α', β'])) # i, α', β'
    L = L_tmp * L_tmp1                                 # α, α', β, β'

    #check whether iMPS is canonical:
    L1 = L * delta(α, α')
    L1mat = array(L1, β, β')
    L1mat0 = Diagonal([L1mat[1,1] for _ in 1:size(L1mat)[1]])

    if norm(L1mat - L1mat0)/norm(L1mat) > ε
        torf = false
        # @show norm(L1mat - L1mat0)
        if showQ
            println("-------------------------------")
            println("Id is NOT L's left eigenvector!")
            @show norm(L1mat - L1mat0)/norm(L1mat)
            # @pt L1mat
        end
    else
        torf = true
        if showQ
            println("⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ⋆ ")
            println("Id is L's left eigenvector!")
            # @pt L1mat
        end
    end
    return R, L, torf
end

function repeatedlyFind(vmat, R, α, β)
    vmat = vmat + vmat'
    vals, vecs, _ = eigsolve(x -> array(ITensor(x, β, β') * R, α, α'), vmat)
    η = vals[1]
    vmat = vecs[1]
    trvmat = tr(vmat)
    vmat *= conj(trvmat)/abs(trvmat)
    max1 = getmaxelm(vmat - vmat')
    return max1, vmat, η
end

function make_canon(Γ, λ, R, L; kwargs...)
    if haskey(kwargs, :cutoffcheck)
        cutoffcheck = kwargs[:cutoffcheck]
    else
        cutoffcheck = 1e-3
    end

    # setup indices
    α = commonind(Γ, Γ, tags="left, bond")
    β = commonind(Γ, Γ, tags="right, bond")

    ## ====================
    ## step 1: find X and Y
    ## ====================
    ### ---------------------
    ### find X
    #### find the dominant right eigenvector (using anonymous function as input)
    ##### note that the eigenvectors returned is already in matrix form!!!
    vals, vecs, _ = eigsolve(
        x -> array(ITensor(x, β, β') * R, α, α'),
        rand(dim(β), dim(β'))
    )
    # if length(vals) > 1 && abs(diff(abs.(vals[1:2]))[1]) < 1e-8
    #     println("degeneracy of VR dominant eigenvalue is found!")
    # end
    η = vals[1]
    η1 = η
    Vᵣmat = vecs[1]
    ##### making eigenvalues of Vᵣ real, if possible
    trVᵣmat = tr(Vᵣmat)
    Vᵣmat *= conj(trVᵣmat)/abs(trVᵣmat)
    max1 = getmaxelm(Vᵣmat - Vᵣmat')
    if max1 < cutoffcheck
        # println("Vᵣmat already is Hermitian!!!!!!!")
    else
        cntmax = 0
        while max1 > cutoffcheck
            cntmax += 1
            max1, Vᵣmat, η = repeatedlyFind(Vᵣmat, R, α, β)

            if max1 < cutoffcheck
                println("Correct Vᵣ found after $cntmax times of repetition")
                break
            end

            if cntmax > 42
                println("NO correct Vᵣ found after $cntmax times of repetition")
                break
            end
        end
   end

    #### switch back to ITensor
    Vᵣ = ITensor(Vᵣmat, β, β') # Vᵣ: β, β'
    if norm(R * Vᵣ - η * replaceinds(Vᵣ, [β, β'] ,[α, α']))/norm(Vᵣ) > cutoffcheck
        println("norm(R * Vᵣ - η * replaceinds(Vᵣ, [β, β'] ,[α, α']))/norm(Vᵣ)
        =
        $(norm(R * Vᵣ - η * replaceinds(Vᵣ, [β, β'] ,[α, α']))/norm(Vᵣ))")
    end

    #### then do the decomposition
    # Vᵣ * W = W * D
    data = eigen(Vᵣmat,sortby=real)
    Dmat = Diagonal(data.values)
    Wmat = data.vectors
    sqrtDmat = sqrt(abs.(Dmat))
    Xmat = Wmat * sqrtDmat
    X = ITensor(Xmat, β, β'')

    ### ---------------------
    ### find Y
    #### find the dominant *left* eigenvector (using anonymous function as input)
    vals, vecs, _ = eigsolve(
        x -> array(ITensor(x, α, α') * L, β, β'),
        rand(dim(α), dim(α'))
    )

    if abs((vals[1] - η)/η) > cutoffcheck
        println("abs(vals[1] - η) = $(abs(vals[1] - η))")
    end

    Vₗmat = vecs[1]
    ##### making eigenvalues of Vₗ real, if possible
    trVₗmat = tr(Vₗmat)
    Vₗmat *= conj(trVₗmat)/abs(trVₗmat)
    max1 = getmaxelm(Vₗmat - Vₗmat')
    if max1 < cutoffcheck
        # println("Vₗmat is already hermitian!!!!!!!!!")
    else
        cntmax = 0
        while max1 > cutoffcheck
            cntmax += 1
            max1, Vₗmat, η = repeatedlyFind(Vₗmat, L, β, α)

            if max1 < cutoffcheck
                println("Correct Vₗ found after $cntmax times of repetition")
                break
            end

            if cntmax > 42
                println("NO correct Vᵣ found after $cntmax times of repetition")
                break
            end
        end
   end
    if abs((η1 - η)/η1) > cutoffcheck
        println("abs((η1 - η)/η1) = $(abs((η1 - η)/η1))")
    end

    #### switch back to ITensor
    Vₗ = ITensor(Vₗmat, α, α')  # Vₗ: α, α'
    if norm(L * Vₗ - η * replaceinds(Vₗ, [α, α'], [β, β']))/norm(Vₗ) > cutoffcheck
        println("norm(L * Vₗ - η * replaceinds(Vₗ, [α, α'], [β, β']))/norm(Vₗ)
        =
        $(norm(L * Vₗ - η * replaceinds(Vₗ, [α, α'], [β, β']))/norm(Vₗ))")
    end

    #### then do the decomposition
    # Vₗ * W = W * D
    data = eigen(Vₗmat,sortby=real)
    Dmat = Diagonal(data.values)
    Wmat = data.vectors
    sqrtDmat = sqrt(abs.(Dmat))
    Xmat = Wmat * sqrtDmat
    Y = ITensor(Xmat, α, α'')

    ## ====================================
    ## step 2: SVD => transpose(Y)λX = Uλ'V to update λ
    ## ====================================
    YtλX = replaceind(Y, α', α) * λ * X # YtλX: uniqidxY, uniqidxX
    U, S, V = svd(YtλX, α''; kwargs...)

    if norm(U * S * V - YtλX)/norm(YtλX) > cutoffcheck
        println("norm(U * S * V - YtλX) / norm(YtλX)
        =
        $(norm(U * S * V - YtλX) / norm(YtλX))")
    end

    αp = commonind(U, S) # U: uniqidxY, αp
    βp = commonind(V, S) # V: βp, uniqidxX
    λ = S
    settags!(λ, "left, bond", αp)
    settags!(λ, "right, bond", βp)
    αp1 = commonind(λ, λ, tags="left, bond")
    βp1 = commonind(λ, λ, tags="right, bond")

    ## ======================================
    ## step 3: construct Γ': V invX Γ invYt U
    ## ======================================
    ### be very careful with the ordering of indices!!!
    Xinv = ITensor(inv(array(X, β, β'')), β'', α)
    Ytinv = ITensor(inv(array(Y, α'', α)), β, α'')
    Γ = V * Xinv * Γ * Ytinv * U
    ### be very careful with the ordering of indices!!!
    replaceinds!(Γ, [αp, βp], [βp1, αp1])

    Γ, λ = normalizeiMPS(Γ, λ, η)
    println("bond dimension =  $(size(λ)[1])")
    return Γ, λ
end
