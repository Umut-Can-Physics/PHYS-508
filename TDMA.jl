# 13/03/2025

# Necessary Modules
using Plots # for plotting
using SpecialFunctions # for erf(x) or error function
using LaTeXStrings # render LaTeX equations on the title of plots

#-#
# FUNCTIONS SECTION
# This section include defined functions to compare solution of differential equation by using TDMA method and exact solution
#-#

"""
Calculate the `A,B,C,D` variables to compute discrete solution of differential problem by using TDMA method.

# Arguments
- `c`,`d`: Start and end points of the 1D grid line 
- `N`: The number of grid points
- `ϵ`: Parameter
"""
function Get_ABCD_Values(c, d, N, ϵ, method)

    # The definition of grid line
    h = (d-c)/(N-1)
    x_range = c:h:d 

    # Initialization of arrays 
    x=zeros(N);a=zeros(N);b=zeros(N);ff=zeros(N)
    for k in 1:N
        x_value = c + (k-1)*h
        x[k] = x_value
        a[k] = x_value
        b[k] = 1
        ff[k] = (1 + ϵ*pi^2)*cos(pi*x[k]) + pi*x[k]*sin(pi*x[k]) 
    end

    # Initialization of arrays 
    A=zeros(N);B=zeros(N);C=zeros(N);F=zeros(N);
    γ=zeros(N);R = (a.*h)./ϵ
    if method==1
        γ = ones(N)
    elseif method==2
        γ = (0.5*abs.(R)).+1
    elseif method==3
        γ = 0.5*R.*coth.(0.5*R)
    end

    for k in 2:N-1
        #A[k] = ϵ/h^2 - a[k]/(2*h)
        A[k] = ϵ/h^2 * (γ[k] - 0.5*R[k])
        B[k] = (2*ϵ*γ[k]) / h^2 + b[k] 
        #C[k] = ϵ/h^2 + a[k]/(2*h)
        C[k] = ϵ/h^2 * (γ[k] + 0.5*R[k])
        F[k] = ff[k]
    end

    # Initial conditions
    B[1] = 1;C[1] = 0;F[1] = -1
    B[N] = 1;A[N] = 0;F[N] = 1

    return A, B, C, F, x, x_range
end

"""
Thomas method to solve differential equations.

# Arguments
- `A`,`B`,`C`,`F`: Variables which are calculated from `Get_ABCD_Values(c, d, N, ϵ)` function defined above
"""
function TDMA(A, B, C, F)

    n = length(F)    

    # Initialization of arrays
    α = zeros(n)
    β = zeros(n)

    α[1] = C[1] / B[1]
    β[1] = F[1] / B[1]

    for i in 2:n-1
        dummy = B[i] - α[i-1] * A[i]
        α[i] = C[i] / dummy
        β[i] = (F[i]+β[i-1]*A[i]) / dummy
    end

    u = zeros(n)

    u[n] = (F[n] + β[n-1] * A[n]) / (B[n]-α[n-1]*A[n])
    for i in n-1:-1:1
        u[i] = α[i] * u[i+1] + β[i]
    end 

    return u
end
"""
Exact solution of the differential equation.

# Arguments
- `x`: 1D coordinates
- `ϵ`: Parameter
"""
function exact(x, ϵ)
    m = x.*erf.(x./sqrt(2*ϵ)).+sqrt(2*ϵ/pi).*exp.(-x.^2 ./(2*ϵ)) 
    n = erf.(1/sqrt(2*ϵ))+sqrt.(2*ϵ/pi)*exp.(-1/(2*ϵ))
    ı = m/n
    return cos.(pi*x) + x + ı
end

#-#
# Results Section
# Thissection includes running of the functions and plot the results
#-#

N=15;ϵ=0.1
c=-1;d=1
method = 1
A, B, C, F, x, x_range = Get_ABCD_Values(c, d, N, ϵ, method)
U = TDMA(A,B,C,F)
plot(x_range, U)
plot!(x_range,exact(x, ϵ),linestyle=:dash)

# Start and end points of grid line
c=-1;d=1

# Effect of the grid refinement
N_list=[10,20,40,80]; ϵ=0.01
P1 = plot()
p1_title = latexstring("\$ ϵ=$(ϵ) \$")
error_list_1 = []
for N in N_list
    A, B, C, F, x, x_range = Get_ABCD_Values(c, d, N, ϵ)
    # Run TDMA method
    U = TDMA(A,B,C,F)
    # Get exact solution
    U_exact = exact(x,ϵ)
    # Error
    errors = push!(error_list_1, abs.(U_exact.-U))
    # ! means hold on the plot to show two of plots simultaneously for comparing exact and discrete solutions
    plot!(P1, x_range,U,title=p1_title,label=latexstring("\$ N=$(N) \$"),xlabel=L"x",dpi=500)
end
savefig(P1, "Effect of grid refinement.png")

# Plot errors for grid refinement
errors_grid_refinement = plot(error_list_1[1], ylabel=L"||e||_\infty", label=latexstring("\$ N=$(N_list[1])\$"),title=L"||e||_\infty=max_{1 \leq i \leq N}|u(x_i)-u_i|", dpi=500)
plot!(error_list_1[2], ylabel=L"||e||_\infty", label=latexstring("\$ N=$(N_list[2])\$"), dpi=500)
plot!(error_list_1[3], ylabel=L"||e||_\infty", label=latexstring("\$ N=$(N_list[3])\$"), dpi=500)
plot!(error_list_1[4], ylabel=L"||e||_\infty", label=latexstring("\$ N=$(N_list[4])\$"), dpi=500)
savefig(errors_grid_refinement,"Errors for Grid Refinement.png")

# Effect of the ϵ refinement
c=-1;d=1
ϵ_list=[1, 1e-1, 1e-2, 1e-3]; N=40
P2 = plot()
p2_title = latexstring("\$ N=$(N) \$")
error_list_2=[]
for ϵ in ϵ_list
    A, B, C, F, x, x_range = Get_ABCD_Values(c, d, N, ϵ)
    # Run TDMA method
    U = TDMA(A,B,C,F)
    # Get exact solution
    U_exact = exact(x,ϵ)
    # Error
    errors = push!(error_list_2, abs.(U_exact.-U))
    # ! means hold on the plot to show two of plots simultaneously for comparing exact and discrete solutions
    plot!(P2, x_range,U,title=p2_title,label=latexstring("\$ ϵ=$(ϵ) \$"),xlabel=L"x",dpi=500)
end
savefig(P2, "Effect of epsilon refinement.png")

# Plot errors for ϵ refinement
errors_epsilon_refinement = plot(error_list_2[1], ylabel=L"||e||_\infty", label=latexstring("\$ ϵ=$(ϵ_list[1])\$"),title=L"||e||_\infty=max_{1 \leq i \leq N}|u(x_i)-u_i|",dpi=500)
plot!(error_list_2[2], ylabel=L"||e||_\infty", label=latexstring("\$ ϵ=$(ϵ_list[2])\$"),dpi=500)
plot!(error_list_2[3], ylabel=L"||e||_\infty", label=latexstring("\$ ϵ=$(ϵ_list[3])\$"),dpi=500)
plot!(error_list_2[4], ylabel=L"||e||_\infty", label=latexstring("\$ ϵ=$(ϵ_list[4])\$"),dpi=500)
savefig(errors_epsilon_refinement,"Errors for Epsilon Refinement.png")

# Get exact solution for a sample parameter
N=41;ϵ=1e-3
_, _, _, _, x, x_range = Get_ABCD_Values(c, d, N, ϵ)
U_exact = exact(x,ϵ)
p_title = latexstring("\$ N=$(N), ϵ=$(ϵ) \$")
# Plot TDMA result and exact solution of the differential equation as a function of x range with the proper plot title
exact_result = plot(x_range, U_exact, linestyle=:dash, label="Exact", linecolor=:black,dpi=500,title=p_title,xlabel=L"x")
savefig(exact_result,"Exact Result.png")