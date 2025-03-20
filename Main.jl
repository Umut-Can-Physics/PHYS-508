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

function input_data(N,c,h)
    ϵ = 0.01
    # Initialization of arrays 
    x=zeros(N);a=zeros(N);b=zeros(N);ff=zeros(N);exact=zeros(N)
    for k in 1:N
        x_value = c + (k-1)*h
        x[k] = x_value
        a[k] = x_value
        b[k] = 1
        ff[k] = (1 + ϵ*pi^2)*cos(pi*x[k]) + pi*x[k]*sin(pi*x[k]) 
        #exact[k] = exact(x,ϵ)
    end
    return ϵ, x, a, b, ff
end

function coeff_matrix(N, ϵ, a, b, ff)
    A=zeros(N);B=zeros(N);C=zeros(N);F=zeros(N);
    for k in 2:N-1
        A[k] = ϵ/h^2 - a[k]/(2*h)
        B[k] = 2*ϵ/h^2 + b[k] 
        C[k] = ϵ/h^2 + a[k]/(2*h)
        F[k] = ff[k]
    end
    # Initial conditions
    B[1] = 1;C[1] = 0;F[1] = -1
    B[N] = 1;A[N] = 0;F[N] = 1
    return A, B, C, F
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