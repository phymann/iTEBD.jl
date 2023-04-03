"""
    initialize iMPS and a tensors
"""
function setupTs(β, J, methodi)
    d = 2
    i = Index(d, "up, site")
    j = Index(d, "left, site")
    k = Index(d, "down, site")
    l = Index(d, "right, site")

    χ = 1
    μ = Index(χ, "left, bond")
    ν = Index(χ, "right, bond")

    #=
        μ-[Γ]-ν
        |
        i
    =#
    Γ = randomITensor(Float64, μ, i, ν)

    #=
        μ-[λ]-ν
    =#
    λ = diagITensor(diag(rand(χ, χ)), μ, ν)

    # make the initial random iMPS canonical
    R, L = canonQ(Γ, λ)
    Γ, λ = make_canon(Γ, λ, R, L)

    Qmat = [
        exp(-β * (-J-2h)) exp(-β * J);
        exp(-β * J) exp(-β * (-J+2h))
    ]

    if methodi == 1
        sqrtQmat = √Qmat
        Qi = ITensor(sqrtQmat, i', i)
        Qj = ITensor(sqrtQmat, j', j)
        Qk = ITensor(sqrtQmat, k, k')
        Ql = ITensor(sqrtQmat, l, l')
    elseif methodi == 2
        # QW = WD
        data = eigen(Qmat)
        sqrtDmat = Diagonal(sqrt.(data.values))
        Wmat = data.vectors
        sqrtQmat = Wmat * sqrtDmat

        Qi = ITensor(sqrtQmat', i', i)
        Qj = ITensor(sqrtQmat', j', j)
        Qk = ITensor(sqrtQmat, k, k')
        Ql = ITensor(sqrtQmat, l, l')
    else
        @error("methodi should be 1 or 2!")
    end

    delta4 = delta(i, j, k, l)

    #=
                i'
                |
                [X]
                [X†]
                |
                i
                |
    j'-[X X†]-j-[a]-l-[X X†]-l'
                |
                k
                |
                [X]
                [X†]
                |
                k'
    =#
    a = Qi * Qj * Qk * Ql * delta4
    replaceinds!(a, [i',j',k',l'], [i, j, k, l])
    return Γ, λ, a
end

"""
    update an iMPS
"""
function updateit(Γ, λ, a)
    μ = commonind(Γ, Γ, tags="left, bond")
    ν = commonind(Γ, Γ, tags="right, bond")
    j = commonind(a, a, tags="left, site")
    l = commonind(a, a, tags="right, site")
    i = commonind(a, a, tags="up, site")
    k = commonind(a, a, tags="down, site")

    it = ITensor(1.)

    #=
    μ-[Γ]-ν
        |
        i
    =#
    it *= Γ

    #=
    μ-[Γ]-ν
        |
    j-[a]-l
        |
        k
    =#
    it *= a

    #=
    μ-[Γa]-νl
        |
        k
    =#
    CL = combiner(j, μ)
    cli = combinedind(CL)
    it *= CL
    settags!(it, "left, bond", cli)

    #=
    μ-[Γa]-ν
        |
        k
    =#
    CR = combiner(l, ν)
    cri = combinedind(CR)
    it *= CR
    settags!(it, "right, bond", cri)

    #=
    μ-[Γa]-ν
        |
        i
    =#
    replaceind!(it, k, i)


    # μj-[λδ]-νl
    λ *= delta(j, l)

    # μ-[λδ]-νl
    λ *= CL
    settags!(λ, "left, bond", cli)
    # μ-[λδ]-ν
    λ *= CR
    settags!(λ, "right, bond", cri)

    return it, λ
end

"""
    get Z from an iMPS
"""
function getZ(Γ, λ, a)
    μ = commonind(Γ, Γ, tags="left, bond")
    ν = commonind(Γ, Γ, tags="right, bond")
    l = commonind(a, a, tags="right, site")
    j = commonind(a, a, tags="left, site")
    i = commonind(a, a, tags="up, site")
    k = commonind(a, a, tags="down, site")
    it1 = ITensor(1.)
    sqrtλ = ITensor(sqrt(array(λ, μ, ν)), μ, ν)
    it1 *= replaceind(sqrtλ, ν, μ') * prime(Γ, "bond") * replaceind(sqrtλ, μ, ν')
    it2 = dag(prime(it1, "bond"))
    replaceind!(it2, i, k)
    vals, _, _, = eigsolve(x -> array(((ITensor(x, ν, l, ν') * it1) * a) * it2, μ, j, μ'), rand(dim(ν), dim(l), dim(ν')))
    return vals[1]
end

"""
    main function for iTEBD
"""
function iTEBDmain(β::Float64, J::Float64, h::Float64; kwargs...)
    # algorithm parameters
    maxdim = get(kwargs, :maxdim, 16)
    cutoff = get(kwargs, :cutoff, 1e-8)
    nrepeat = get(kwargs, :nrepeat, 420)
    cutoffcheck = get(kwargs, :cutoffcheck, 1e-2)
    methodi = get(kwargs, :methodi, 2)
    ## method for obtain a tensor
    #=
        1: simply W = √Q, W W = Q
        2: first eigen decomposition, then obtain W W† = Q
    =#

    # prepare tensors
    # ---------------------------------------
    Γ, λ, a = setupTs(β, J, methodi)

    # update tensors
    # ---------------------------------------
    for counti in 1:nrepeat
        @show counti
        Γold = Γ
        λold = λ

        Γ, λ = updateit(Γold, λold, a)

        R, L = canonQ(Γ, λ)
        Γ, λ = make_canon(Γ, λ, R, L;
        maxdim = maxdim,
        ishermitian = true,
        cutoffcheck = cutoffcheck,
        cutoff = cutoff)

        if size(λ) == size(λold)
            cond = norm(matrix(λ) - matrix(λold))
            println("convergence criterion = $cond")
            if cond < cutoff
                break
            end
        end
    end

    if abs(h) < 1e-5
        FE = -log(getZ(Γ, λ, a))/β
        FEe = ising_free_energy(β, J)
        println("------------------------------")
        println("relative free energy diff:")
        println("$(abs((FE - FEe)/FEe))")
    end
    return FE
end
