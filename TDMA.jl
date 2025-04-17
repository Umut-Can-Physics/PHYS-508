# Necessary Modules
using Plots # for plotting
using SpecialFunctions # for erf(x) or error function
using LaTeXStrings # render LaTeX equations on the title of plots
using Revise
includet("Main.jl")

#-#
# Results Section
# This section includes running of the functions and plot the results
#-#

N=108;ϵ=0.01
c=-1;d=1
h = (d-c)/(N-1)
x_range = c:h:d 

param = "linear"
u = 0
method = 1
A, B, C, F, x, x_range, u_exact_list_1 = CoeffMatrix(c, h, N, ϵ, u, param, method)
U1 = TDMA(A,B,C,F)

method = 2
A, B, C, F, x, x_range, u_exact_list_2 = CoeffMatrix(c, h, N, ϵ, u, param, method)
U2 = TDMA(A,B,C,F)

method = 3
A, B, C, F, x, x_range, u_exact_list_3 = CoeffMatrix(c, h, N, ϵ, u, param, method)
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
Max_errors_plot = plot(N_list,RMS[1,:],xlabel=L"N",ylabel=L"||e||_\infty", title=plt_titlee, label="Centeral Difference", dpi=500,
markershape=:circle)
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

