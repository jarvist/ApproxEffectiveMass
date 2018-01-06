using ApproxFun
println("Hello world! Imports successful...")


#f=Fun(sin,Fourier(-2*pi..2*pi)) # Test Approxfun call.

function FunEffectiveMass(f)
    extrema=roots(f')

    println("Extrema look like: {k_x}=", extrema)
    println("Energy at Extrema: f({k_x})=",f(extrema) )
    println("f''({k_x}) at Extrema: =",f''(extrema) )

    # NOT FULLY WORKING YET; BUT GETTING THERE!
    a=5.431E-10 # 5.431 Angstrom
    kmax=10*pi/a # Factor of 10 is magic number for discretisation (# of k-points) in read in EIGENVAL
    hbar = 1.0545718E-34 # hbar; SI 
    me = 9.10938356E-31 # Mass of Electron, SI, kg
    q=1.602E-19 # Electron charge; to covert from VASP eV --> SI, Joules
    
    meff=Float64[] # A bit inelegant; certainly a more sophisticated way to do this in Julia
    for d2Edk2 in f''(extrema)
        push!(meff,  q * kmax^2 * hbar / ( d2Edk2 * me)) # Not sure of this; many text books drop odd factors :^)
    end 

    println("Effective masses: ",meff) 

    plot!(f)
    scatter!(extrema,f(extrema);color=:green)
    ytextoffset=-0.8
    for ex in extrema
        annotate!(ex,f(ex)+ytextoffset,text(@sprintf("%.2f",f''(ex)),8,:center,:purple)) # Size 10 ?
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

    plot(bands,title="Raw Band structure; from EIGENVAL",xaxis="Forever lost in the Brillouin Zone", yaxis="Dos (eV)")

    return bands
end

# From: https://github.com/ApproxFun/ApproxFun.jl/issues/275 , courtesy of private communication with Sheehan Olver
# Least squares approximation of data on an evenly spaced grid with Chebyshev series
function vandermonde(S,n,x::AbstractVector)
    V=Array(Float64,length(x),n)
    for k=1:n
        V[:,k]=Fun(S,[zeros(k-1);1])(x)
    end
    V
end

# For ...(this)... case, make sure `length(pts) >> n`.
function ApproxFunVandermonde(vals,n=20,  lower=0.0, upper=360.0)
    c=Fourier(lower..upper) #Define Fourier domain in this range (to match data imported)

    pts=collect(1.0:length(vals)) # Points 1..length(vals)

    V=vandermonde(c,n,pts)
#    println(V,pts,vals)
#    print(V\vals)
    # Are you ready for the magic?
    af=Fun(c,V\vals) # Approximate Function (af)
    # me is now an ApproxFun representation of the tabulated data.
    # As a Chebyshev polynomial fit we can do all sorts of differentiation + integration.
    return af
end

# Functions defined.
# Here is the more 'script' part of the code which actually runs.

using Plots
bands=read_EIGENVAL(open("EIGENVAL","r"))

plt = plot(title = "ApproxFun Effective Masses", xaxis = "Forever lost in the Brillouin Zone", yaxis="DoS (eV)")

for band in bands
    ApproxFunBand=ApproxFunVandermonde(band,20,1,40)
    println("\nBand...",band)
    FunEffectiveMass(ApproxFunBand)
end

png("plot.png") # Save last plot produced...

