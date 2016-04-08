using ApproxFun
println("Hello world! Imports successful...")


f=Fun(sin,Fourier([-2*pi,2*pi]))

extrema=roots(f')

println("Extrema look like: (x)", extrema)
println("f'' at Extrema: ",f''(extrema) )

using Plots
#plot(f)
#scatter!(extrema,f(extrema);color=:green)
#show()

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

    band_energies = [[] for i in 1:nbands] # Initialise a load of empty set, one per band

    for i in 1:npoints
        x,y,z,weight = split(readline(f))
        for j in 1:nbands
            fields = split(readline(f))
            id = parse(Int,fields[1])
            energy = parse(Float64,fields[2])
            push!(band_energies[id],energy)
        end
        empty=readline(f)
    end

    println(band_energies[1])

    plot(band_energies)
end

read_EIGENVAL(open("EIGENVAL","r"))

show() # Show any plots produced...
