# Solving u_xx = -6x
# u(0) = 0, u(1) = 1

using Plots
using LinearAlgebra

ys = [2^2, 2^3, 2^4, 2^5, 2^6]

hs = 1 ./ys

i = 4
h = hs[i]

x = range(0,step=h,stop=1)

N = length(x)-2



# We use hat functions as base functions
# φᵢ(x) = 1/h (x - (i-1)h) ∀ (i-1)h ≤ x ≤ ih
# φᵢ(x) = -1/h (x - (i+1)h) ∀ ih ≤ x ≤ (i+1)h

# Assembly row-wise
A = zeros(N,N)

for i in range(1,stop=N)
    A[i,i] = 2/h;
end

for i in range(1,stop=N-1)
    A[i,i+1] = -1/h;
end

for i in range(2,stop=N)
    A[i,i-1] = -1/h;
end


# Assembly element-wise
A = zeros(N,N)

Ae = [1/h -1/h; -1/h 1/h]

for i in 2:N
    A[i-1:i,i-1:i] = A[i-1:i,i-1:i] .+ Ae
end

A[1,1] = A[1,1] + 1/h
A[N,N] = A[N,N] + 1/h

# b = 12h*ones(N)

b = zeros(N)
for i in range(1,N)
    b[i] = h^2*(6*i-4)   # analytical integration for b yields this result
end

num_sol = A\b
plot(x,vcat(0,num_sol,0))


analy_sol = -x .^3 .+x  # analytical solution
plot(analy_sol)





# A = zeros(N,N)
#
# Ae = [1/h -1/h; -1/h 1/h]
#
# for i in 2:N
#     A[i-1:i,i-1:i] = A[i-1:i,i-1:i] .+ Ae
# end
#
# A[1,1] = A[1,1] + 1/h
# A[N,N] = A[N,N] + 1/h
#
# num_sol_2 = A\b


function convergence_test(n) # n is the number of different meshes
    errs = zeros(n)
    for i in range(1,n)
        h = hs[i]

        x = range(0,step=h,stop=1)

        N = length(x)-2

        h = hs[i]

        x = range(0,step=h,stop=1)

        N = length(x)-2 # Number of finite elements

        # Use hat functions


        A = zeros(N,N)

        Ae = [1/h -1/h; -1/h 1/h]

        for i in 2:N
            A[i-1:i,i-1:i] = A[i-1:i,i-1:i] .+ Ae
        end

        A[1,1] = A[1,1] + 1/h
        A[N,N] = A[N,N] + 1/h

        # b = 12h*ones(N)

        b = zeros(N)
        for i in range(1,stop=N)
            b[i] = h^2*(6*i-4)
        end

        num_sol = A\b
        num_sol = vcat(0,num_sol,0)
        analy_sol = -x .^3 .+x

        err = h*sqrt(norm(num_sol-analy_sol))
        errs[i] = log(2,err)
    end
    return errs
end

n = 5
errs = convergence_test(n)
plot(range(1,stop=n),errs)   # convergence_test shows accuracy order = 1
