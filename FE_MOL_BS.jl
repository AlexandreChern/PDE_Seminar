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
for j in 1:N-1
	A[j,j] = 1
end

for p in 1:N-2
	A[p+1,p] = -1
end

A = A*(-a/h)

b = zeros(N-1,1)


u = sin.(pi*x[2:N])

w = [c*sin(t[1]); u]
plot(x,w)
pause(1)
clf()

for i in 2:M

	b[1] = c*(a/h)*sin(t[i-1])
	 
	u[:] = u[:] .+ k.*A*u[:] .+ k*b
	w = [c*sin(t[i]); u]

	exact = zeros(Nfine)
	for v in 1:Nfine
		if xfine[v] - a*t[i] > 0
			exact[v] = sin.(pi*(xfine[v]-a*t[i]))
		else
			exact[v] = c*sin(t[i] - xfine[v]/a)
		end
	end
	plot(x,w,xfine,exact)
	pause(1)
	clf()

end




