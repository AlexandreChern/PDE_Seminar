using Plots
#using PyPlot

h = 0.01
lambda = .5
k = lambda*h

T = 2
a = 0.5
c = 0

x = 0:h:1
t = 0:k:T


xfine = 0:h/20:1
Nfine = length(xfine)


N = length(x)
M = length(t)

A = zeros(N-1,N-1)
A1 = zeros(N-1,N-1)
for j in 1:N-1
    A1[j,j] = 1
end

A2 = zeros(N-1,N-1)
for j in 1:N-1
    A2[j,j] = 3
end
for j in 1:N-2
    A2[j,j+1] = 1
end

for j in 2:N-1
    A2[j,j-1] = -4
end

A3 = zeros(N-1,N-1)

for j in 1:N-1
    A3[j,j] = 1
end
for j in 1:N-2
    A3[j,j] = 1
end
for j in 2:N-1
    A3[j,j-1] = -2
end

A = A1 + (-a*k/(2h))*A2 + (a^2*k^2/(2*h^2))*A3

b = zeros(N-1,1)

u = sin.(pi*x[2:N])

w = [c*sin.(t[1]);u]

plot(x,w)


for i in 2:M
    b[1] = c*(a/h)*sin(t[i-1])
    u[:] = u[:] .+ k.*A*u[:] .+ k*b
    w = [c*sin(t[i]);u]

    exact = zeros(Nfine)
    for v in 1:Nfine
        if xfine[v] - a*t[i] > 0
            exact[v] = sin.(pi*(xfine[v] - a*t[i]))
        else
            exact[v] = c*sin(t[i] - xfine[v]/a)
        end
    end
    plot(x,w)
    pause(1)
    clf()
end
