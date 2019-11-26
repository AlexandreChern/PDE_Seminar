#Julia script for solving advection u_t + au_x = 0 with a > 0, f(x) = sin(pi*x) and g(t) = c*sin(t)
#Applying backward difference in space (applying MOL), forward-Euler in time on domain x in [0, 1]
#and temporal domain t in [0, T]

#to call this script from the Julia repl, type include("FE_MOL_BS.jl")

# h:      grid spacing
# lambda: h/h
# k:      time step
# xfine:  fine grid on which to plot exact solution



using Pkg
#Pkg.add("PyPlot")
using PyPlot
using LinearAlgebra

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

D = zeros(N,N)

for i in 2:N-1
    D[i,i+1] = 1/(2*h)
    D[i,i-1] = -1/(2*h)
end

D[1,1] = -1/h
D[1,2] = 1/h

D[N,N-1] = -1/h
D[N,N] = 1/h

Hinv = h*Matrix{Float64}(I,N,N)
Hinv[1,1] = 2/h
Hinv[N,N] = 2/h






b = zeros(N-1,1)
sig = -a/20

u = sin.(pi*x[2:N])
e1 = zeros(N,1)
e1[1] = 1
plot(x,u)
pause(1)
clf()

g(t) = c*sin(t)

for i in 2:M	 
	u[:] = u[:] .- a.*k.*D*u[:] .+ sig*Hinv*(u[1] - g(t[i-1]))*e1

	exact = zeros(Nfine)
	for v in 1:Nfine
		if xfine[v] - a*t[i] > 0
			exact[v] = sin.(pi*(xfine[v]-a*t[i]))
		else
			exact[v] = c*sin(t[i] - xfine[v]/a)
		end
	end
	plot(x,u,xfine,exact)
	pause(1)
	clf()

end




