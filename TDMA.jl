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
function Get_ABCD_Values(c, h, N, ϵ, method)

    # Initialization of arrays 
    x=zeros(N);a=zeros(N);b=zeros(N);ff=zeros(N)
    u_exact_list=[]
    for k in 1:N
        x_value = c + (k-1)*h
        x[k] = x_value
        a[k] = x_value
        b[k] = 1
        push!(u_exact_list, exact(x_value, ϵ))
        ff[k] = (1 + ϵ*pi^2)*cos(pi*x[k]) + pi*x[k]*sin(pi*x[k]) 
    end

    # Initialization of arrays 
    A=zeros(N);B=zeros(N);C=zeros(N);F=zeros(N);
    γ=zeros(N);R=(a.*h)./ϵ
    if method==1
        γ = ones(N)
    elseif method==2
        γ = (0.5*abs.(R)).+1
    elseif method==3
        for k in 2:N-1
            if R[k]>1e-6
                γ[k] = 0.5*R[k].*coth.(0.5*R[k])    
            else
                γ[k] = 1
            end
        end
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

    return A, B, C, F, x, x_range, u_exact_list
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

N=8;ϵ=0.01
c=-1;d=1
h = (d-c)/(N-1)
x_range = c:h:d 

method = 1
A, B, C, F, x, x_range, u_exact_list_1 = Get_ABCD_Values(c, h, N, ϵ, method)
U1 = TDMA(A,B,C,F)

method = 2
A, B, C, F, x, x_range, u_exact_list_2 = Get_ABCD_Values(c, h, N, ϵ, method)
U2 = TDMA(A,B,C,F)

method = 3
A, B, C, F, x, x_range, u_exact_list_3 = Get_ABCD_Values(c, h, N, ϵ, method)
U3 = TDMA(A,B,C,F)
plt_title=latexstring("\$ N=$(N), ϵ=$(ϵ) \$")
Solutions = plot(x_range, [U1, U2, U3], xlabel=L"x", title=plt_title, dpi=500)

# Fine grid for exact solution
N_fine=400
h_fine = (d-c)/(N_fine-1)
x_range_fine = c:h_fine:d 
x_value = [c + (k-1)*h_fine for k in 1:N_fine]
plot(x_range_fine, exact(x_value, ϵ), linestyle=:dash, label="Exact", dpi=500)
savefig(Solutions,"Solutions for N=$(N).png")


# PLOT ROUND MEAN SQUARED ERRORS (GRID REFINEMENT)
ϵ=0.01
RMS = zeros(3,4)
N_list = [10,20,40,80]
for (idx,N) in enumerate(N_list)

    c=-1;d=1
    h = (d-c)/(N-1)
    x_range = c:h:d 

    method = 1
    A, B, C, F, x, x_range, u_exact_list_1 = Get_ABCD_Values(c, h, N, ϵ, method)
    U1 = TDMA(A,B,C,F)
    e2_1 = sqrt((sum((U1.-u_exact_list_1).^2)) / N) 
    RMS[1,idx] = e2_1

    method = 2
    A, B, C, F, x, x_range, u_exact_list_2 = Get_ABCD_Values(c, h, N, ϵ, method)
    U2 = TDMA(A,B,C,F)
    e2_2 = sqrt((sum((U2.-u_exact_list_2).^2)) / N)
    RMS[2,idx] = e2_2

    method = 3
    A, B, C, F, x, x_range, u_exact_list_3 = Get_ABCD_Values(c, h, N, ϵ, method)
    U3 = TDMA(A,B,C,F)
    e2_3 = sqrt(sum((U3.-u_exact_list_3).^2) / N) 
    RMS[3,idx] = e2_3
end
plt_titlee = latexstring("\$ ϵ=$(ϵ)\$")
RMS_plot = plot(N_list,RMS[1,:],xlabel=L"N",ylabel=L"||e||_2", title=plt_titlee, label="Centeral Difference", dpi=500)
plot!(N_list,RMS[2,:],label="Directed Difference",dpi=500)
plot!(N_list,RMS[3,:],label="Exponential",dpi=500)
savefig(RMS_plot,"RMS Plot.png")

# PLOT MAXIMUM ERRORS (GRID REFINEMENT)
ϵ=0.01
Max_errors = zeros(3,4)
N_list = [10,20,40,80]
for (idx,N) in enumerate(N_list)

    c=-1;d=1
    h = (d-c)/(N-1)
    x_range = c:h:d 

    method = 1
    A, B, C, F, x, x_range, u_exact_list_1 = Get_ABCD_Values(c, h, N, ϵ, method)
    U1 = TDMA(A,B,C,F)
    e_max_1 = maximum(abs.(U1.-u_exact_list_1))
    RMS[1,idx] = e_max_1

    method = 2
    A, B, C, F, x, x_range, u_exact_list_2 = Get_ABCD_Values(c, h, N, ϵ, method)
    U2 = TDMA(A,B,C,F)
    e_max_2 = maximum(abs.(U2.-u_exact_list_2))
    RMS[2,idx] = e_max_2

    method = 3
    A, B, C, F, x, x_range, u_exact_list_3 = Get_ABCD_Values(c, h, N, ϵ, method)
    U3 = TDMA(A,B,C,F)
    e_max_3 = maximum(abs.(U3.-u_exact_list_3))
    RMS[3,idx] = e_max_3
end
plt_titlee = latexstring("\$ ϵ=$(ϵ)\$")
Max_errors_plot = plot(N_list,RMS[1,:],xlabel=L"N",ylabel=L"||e||_\infty", title=plt_titlee, label="Centeral Difference", dpi=500)
plot!(N_list,RMS[2,:],label="Directed Difference",dpi=500)
plot!(N_list,RMS[3,:],label="Exponential",dpi=500)
savefig(Max_errors_plot,"Max Errors Plot.png")

# PLOT ROUND MEAN SQUARED ERRORS (EPSILON REFINEMENT)
RMS_epsilon = zeros(3,4)
ϵ_list = [1,0.1,0.01,0.001]
N = 40
for (idx,ϵ) in enumerate(ϵ_list)

    c=-1;d=1
    h = (d-c)/(N-1)
    x_range = c:h:d 

    method = 1
    A, B, C, F, x, x_range, u_exact_list_1 = Get_ABCD_Values(c, h, N, ϵ, method)
    U1 = TDMA(A,B,C,F)
    e2_1 = sqrt((sum((U1.-u_exact_list_1).^2)) / N) 
    RMS_epsilon[1,idx] = e2_1

    method = 2
    A, B, C, F, x, x_range, u_exact_list_2 = Get_ABCD_Values(c, h, N, ϵ, method)
    U2 = TDMA(A,B,C,F)
    e2_2 = sqrt((sum((U2.-u_exact_list_2).^2)) / N)
    RMS_epsilon[2,idx] = e2_2

    method = 3
    A, B, C, F, x, x_range, u_exact_list_3 = Get_ABCD_Values(c, h, N, ϵ, method)
    U3 = TDMA(A,B,C,F)
    e2_3 = sqrt(sum((U3.-u_exact_list_3).^2) / N) 
    RMS_epsilon[3,idx] = e2_3
end
plt_titlee = latexstring("\$ N=$(N)\$")
RMS_epsilon_plot = plot(ϵ_list, RMS_epsilon[1,:], xlabel="ϵ(log)", ylabel=L"||e||_2", title=plt_titlee, 
label="Centeral Difference", dpi=500,
xticks=ϵ_list,xaxis=:log)
plot!(ϵ_list,RMS_epsilon[2,:],label="Directed Difference",dpi=500)
plot!(ϵ_list,RMS_epsilon[3,:],label="Exponential",dpi=500)
savefig(RMS_epsilon_plot,"RMS Epsilon Plot.png")

# PLOT MAXIMUM ERRORS (Epsilon REFINEMENT)
Max_errors_epsilon = zeros(3,4)
ϵ_list = [1,0.1,0.01,0.001]
N = 40
for (idx,ϵ) in enumerate(ϵ_list)

    c=-1;d=1
    h = (d-c)/(N-1)
    x_range = c:h:d 

    method = 1
    A, B, C, F, x, x_range, u_exact_list_1 = Get_ABCD_Values(c, h, N, ϵ, method)
    U1 = TDMA(A,B,C,F)
    e_max_1 = maximum(abs.(U1.-u_exact_list_1))
    Max_errors_epsilon[1,idx] = e_max_1

    method = 2
    A, B, C, F, x, x_range, u_exact_list_2 = Get_ABCD_Values(c, h, N, ϵ, method)
    U2 = TDMA(A,B,C,F)
    e_max_2 = maximum(abs.(U2.-u_exact_list_2))
    Max_errors_epsilon[2,idx] = e_max_2

    method = 3
    A, B, C, F, x, x_range, u_exact_list_3 = Get_ABCD_Values(c, h, N, ϵ, method)
    U3 = TDMA(A,B,C,F)
    e_max_3 = maximum(abs.(U3.-u_exact_list_3))
    Max_errors_epsilon[3,idx] = e_max_3
end
plt_titlee = latexstring("\$ N=$(N)\$")
Max_errors_epsilon_plot = plot(ϵ_list, Max_errors_epsilon[1,:],xlabel="ϵ(log)",ylabel=L"||e||_\infty", title=plt_titlee, 
label="Centeral Difference", dpi=500,
xticks=ϵ_list,xaxis=:log)
plot!(ϵ_list,Max_errors_epsilon[2,:],label="Directed Difference",dpi=500)
plot!(ϵ_list,Max_errors_epsilon[3,:],label="Exponential",dpi=500)
savefig(Max_errors_epsilon_plot,"Max Errors Epsilon Plot.png")