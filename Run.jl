using Revise
using SpecialFunctions
using Plots

includet("Main.jl")

# The definition of grid line
N=41
c=-1;d=1
h = (d-c)/(N-1)
x_range = c:h:d 

# Functions
ϵ, x, a, b, ff = input_data(N, c, h)
A, B, C, F = coeff_matrix(N, ϵ, a, b, ff)
U = TDMA(A, B, C, F)

plot(U, x_range)
plot(u_exact, x_range)