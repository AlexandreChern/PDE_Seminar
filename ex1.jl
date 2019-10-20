# Solve advection with initial data f(x) = sin(pi*x) and g(t) =
h = 0.1

lambda = 0.5

k = h*lambda
a = 1

T = 2

x = 0:h:1
t = 0:k:T

N = length(x)
M = length(t)


A = zeros(N-1,N-1)

for j in 1:N-1
    A[j,j] = 1/h
end

for j in 2:N-1
    A[j,j-1] = -1/h
end

b = zeros(N-1)

u = sin.(pi*x[2:N])

for j in 2:M
    u[:] = u[:] - a*k*A*u[:] + k*b
end

for j in 2:M
    u .= u .- a*k*A*u + k*b
end


for j in 2:M
    u = u - a*k*A*u + k*b ## This doesn't work
end
