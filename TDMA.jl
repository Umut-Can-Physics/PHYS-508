# Page 42, Ex. 4.4

using Plots
using SpecialFunctions

N = 101
c = -1; d = 1
h = (d-c)/(N-1)
x_range = c:h:d 

ϵ = 10
x=zeros(N);a=zeros(N);b=zeros(N);ff=zeros(N)
for k in 1:N
    x_value = c + (k-1)*h
    x[k] = x_value
    a[k] = x_value
    b[k] = 1
    ff[k] = (1 + ϵ*pi^2)*cos(pi*x[k]) + pi*x[k]*sin(pi*x[k]) 
end

A=zeros(N);B=zeros(N);C=zeros(N);F=zeros(N);
for k in 2:N-1
    A[k] = -ϵ/h^2 - a[k]/(2*h)
    B[k] = 2*ϵ/h^2 + b[k] 
    C[k] = -ϵ/h^2 + a[k]/(2*h)
    F[k] = ff[k]
end

B[1] = 1
C[1] = 0
F[1] = -1

B[N] = 1
A[N] = 0
F[N] = 1

# Thomas method
function TDMA(A, B, C, F)
    n = length(F)    

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

U = TDMA(A,B,C,F)

plot(x_range,U)

function exact(x,ϵ)
    m = x.*erf.(x./sqrt(2*ϵ)).+sqrt(2*ϵ/pi).*exp.(-x.^2 ./(2*ϵ)) 
    n = erf.(1/sqrt(2*ϵ))+sqrt.(2*ϵ/pi)*exp.(-1/(2*ϵ))
    ı = m/n
    return cos.(pi*x) + x + ı
end
U_exact = exact(x,ϵ)
plot!(x_range,U_exact, linestyle=:dash)