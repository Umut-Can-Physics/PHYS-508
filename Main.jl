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