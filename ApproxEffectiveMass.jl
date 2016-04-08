using ApproxFun

println("Hello world! Imports successful...")

f=Fun(sin,Fourier([-2*pi,2*pi]))

extrema=roots(f')

println("Extrema look like: (x)", extrema)
println("f'' at Extrema: ",f''(extrema) )


using Plots
plot(f)
scatter!(extrema,f(extrema);color=:green)
show()
