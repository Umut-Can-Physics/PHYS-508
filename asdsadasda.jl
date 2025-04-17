using Plots
using SpecialFunctions
using LaTeXStrings

"""
Compute coefficients for TDMA discretization of linear or nonlinear BVP.
"""
function CoeffMatrix(c, h, N, ϵ, u, param, method)
    # Spatial grid
    x_range = range(c, stop = c + (N-1)*h, length = N)
    x = collect(x_range)
    # Preallocate tridiagonal arrays
    A = zeros(N); B = zeros(N); C = zeros(N); F = zeros(N)

    if param == "linear"
        error("Linear branch not implemented in this snippet.")
    elseif param == "non-linear"
                # Prepare nonlinear terms
        # a(x) same advection coefficient as in linear case
        a = copy(x)
        b = zeros(N)
        # Compute source and Jacobian from current iterate u
        for k in 2:N-1
            g = u[k] + u[k]^2             # g(u) = u + u^2
            dg = 1 + 2*u[k]               # dg/du
            b[k] = dg                    
            f0 = exp(-2 * x[k] / sqrt(ϵ)) # RHS term
            F[k] = f0 - g + dg*u[k]
        end
        # Scharfetter–Gummel weights for advection (using a = x)
        R = a .* h ./ ϵ
        gamma = if method == 1
            ones(N)
        elseif method == 2
            0.5 * abs.(R) .+ 1
        else
            tmp = ones(N)
            for k in 2:N-1
                if abs(R[k]) > 1e-6
                    tmp[k] = 0.5*R[k]*coth(0.5*R[k])
                end
            end
            tmp
        end
        # Build coefficients
        for k in 2:N-1
            A[k] = ϵ/h^2*(gamma[k] - 0.5*R[k])
            B[k] = 2*ϵ*gamma[k]/h^2 + b[k]
            C[k] = ϵ/h^2*(gamma[k] + 0.5*R[k])
        end
        # Boundary conditions
        B[1], C[1], F[1] = 1, 0, u[1]
        B[N], A[N], F[N] = 1, 0, u[N]
    else
        error("Unknown param: $param")
    end

    return A, B, C, F, x, x_range
end

"""
Thomas algorithm (TDMA) for tridiagonal systems.
"""
function TDMA(A, B, C, F)
    n = length(F)
    α = zeros(n); β = zeros(n)
    α[1] = C[1]/B[1]; β[1] = F[1]/B[1]
    for i in 2:n-1
        denom = B[i] - α[i-1]*A[i]
        α[i] = C[i]/denom
        β[i] = (F[i] + β[i-1]*A[i]) / denom
    end
    u = zeros(n)
    u[n] = (F[n] + β[n-1]*A[n]) / (B[n] - α[n-1]*A[n])
    for i in n-1:-1:1
        u[i] = α[i]*u[i+1] + β[i]
    end
    return u
end

"""
Initial guess for nonlinear solver.
"""
function UInitial(x, ϵ)
    return 1 .+ (exp(-1/sqrt(ϵ)) - 1) .* x
end

"""
Solve nonlinear BVP by fixed-point/Newton linearization and TDMA.
Returns solution vector and corresponding x-range.
"""
function solve_nonlinear(c, d, N, ϵ, method; tol=1e-6, maxiter=100)
    h = (d - c) / (N - 1)
    x_range = range(c, stop = d, length = N)
    u_prev = UInitial(collect(x_range), ϵ)
    for iter in 1:maxiter
        A, B, C, F, x_vals, _ = CoeffMatrix(c, h, N, ϵ, u_prev, "non-linear", method)
        u_new = TDMA(A, B, C, F)
        err_rel = maximum(abs.(u_new .- u_prev)) / maximum(abs.(u_prev))
        @info "Iteration $iter: err_rel = $err_rel"
        u_prev .= u_new
        if err_rel < tol
            break
        end
    end
    return u_prev, x_range
end

N = 108; c = 0; d = 1; ϵ = 1; method = 2
sol, xs = solve_nonlinear(c, d, N, ϵ, method)
plot(xs, sol, label = "Numerical")
plot!(xs, exp.(-xs ./ sqrt(ϵ)), label = "Exact")