# Necessary Modules
using Plots # for plotting
using SpecialFunctions # for erf(x) or error function
using LaTeXStrings # render LaTeX equations on the title of plots
using Revise
includet("Main.jl")

N=108
c=0;d=1
h = (d-c)/(N-1)
x_range = c:h:d 
x_value = [c + (k-1)*h for k in 1:N]
ϵ = 1

ExactSol = Exact_2(x_value, ϵ)
u_prev = UInitial(x_value, ϵ)

method = 2
param = "non-linear"

tol=1e-6; maxiter=100;
for iter in 1:maxiter
    A, B, C, F, x, x_range, _ = CoeffMatrix(c, h, N, ϵ, u_prev, param, method)
    u_new = TDMA(A,B,C,F)
    err_rel = maximum(abs.(u_new .- u_prev)) / maximum(abs.(u_prev))
    @info "Iteration $iter: err_rel = $err_rel"
    u_prev .= u_new
    if err_rel < tol
        break
    end
end

plot(x_range, u_prev)
plot!(x_range, ExactSol)