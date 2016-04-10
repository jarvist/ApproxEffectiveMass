using ApproxFun
using Plots
println("Hello world! Imports successful...")


f=Fun(sin,Fourier([-2*pi,2*pi]))

function FunEffectiveMass(f)
    extrema=roots(f')

    println("Extrema look like: (x)", extrema)
    println("f'' at Extrema: ",f''(extrema) )

    plot!(f)
    scatter!(extrema,f(extrema);color=:green)
    ytextoffset=0.1
    for ex in extrema
        annotate!(ex,f(ex)+ytextoffset,text(@sprintf("%.2f",f''(ex)),:center,:orange))
    end
end

#FunEffectiveMass(f)

# Code to read in horrid VASP EIGENVAL form, with understanding file format from reading jkitchin's JASP:
# https://github.com/jkitchin/jasp/blob/31eda6cc3e64d6e953d9d372b80abb2d3e76559c/jasp/jasp_bandstructure.py#L66-L107

# From my SMASH.jl
function readnlines(f,n)
    local lines=""
    local i=1
    for i=1:n
        lines=lines*readline(f)
    end
    return (lines)
end

#readmatrix(f, nlines) = readdlm(IOBuffer(string([readline(f) for i in 1:nlines])))
readmatrix(f, nlines) = readdlm(IOBuffer(readnlines(f,nlines)))

function read_EIGENVAL(f::IOStream)
    junk=readnlines(f,5)
    unknown, npoints, nbands = split(readline(f))
    npoints=parse(Int,npoints)
    nbands=parse(Int,nbands)
    println("unknown: $(unknown) npoints: $(npoints) nbands: $(nbands) ")
    empty=readline(f)

    bands = [Float64[] for i in 1:nbands] # Initialise a load of empty set, one per band

    for i in 1:npoints
        x,y,z,weight = split(readline(f))
        for j in 1:nbands
            fields = split(readline(f))
            id = parse(Int,fields[1])
            energy = parse(Float64,fields[2])
            
            push!(bands[id],energy)
        end
        empty=readline(f)
    end

    println(bands[1])

    plot(bands)

    return bands
end

# From: https://github.com/ApproxFun/ApproxFun.jl/issues/275 , courtesy of private communication with Sheehan Olver
# Least squares approximation of data on an evenly spaced grid with Chebyshev series
function vandermonde(S,n,x::AbstractVector)
    V=Array(Float64,length(x),n)
    for k=1:n
        V[:,k]=Fun([zeros(k-1);1],S)(x)
    end
    V
end

# For ...(this)... case, make sure `length(pts) >> n`.
function ApproxFunVandermonde(vals,n=20,  lower=0.0, upper=360.0)
    c=Fourier([lower,upper]) #Define Fourier domain in this range (to match data imported)

    pts=collect(1.0:length(vals)) # Points 1..length(vals)

    V=vandermonde(c,n,pts)
    println(V,pts,vals)
    print(V\vals)
    # Are you ready for the magic?
    af=Fun(V\vals,c) # Approximate Function (af)
    # me is now an ApproxFun representation of the tabulated data.
    # As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.
    return af
end

bands=read_EIGENVAL(open("EIGENVAL","r"))

plot()
for band in bands
    ApproxFunBand=ApproxFunVandermonde(band,20,1,40)
    FunEffectiveMass(ApproxFunBand)
end

show() # Show any plots produced...
