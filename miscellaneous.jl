function jwdisplay(x,str)
    println("---------------------------")
    println("$str is ")
    display(x)
end

function jwchk(x::Bool)
    if !x
        error("\n
        ================================
        sanity check failed
        ================================\n")
    end
end

function ising_free_energy(β::Real, J::Real=1.0)
    k = β * J
    c = cosh(2 * k)
    s = sinh(2 * k)
    xmin = 0.0
    xmax = π
    integrand(x) = log(c^2 + √(s^4 + 1 - 2 * s^2 * cos(2x)))
    integral, err = quadgk(integrand, xmin, xmax)::Tuple{Float64,Float64}
    return -(log(2) + integral / π) / (2 * β)
end

function getmaxelm(mat)
    return maximum(abs.(mat))
end
