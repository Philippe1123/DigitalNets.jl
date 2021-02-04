


struct DigitalNet32{s} <: AbstractDigitalNets{s}
    C::Matrix{UInt32} # generating matrix
    n::Int64 # max number of points in the lattice rule
    recipid::Float64 # internal

end

struct DigitalNet64{s} <: AbstractDigitalNets{s}
    C::Matrix{UInt64} # generating matrix
    n::Int64 # max number of points in the lattice rule
    recipid::Float64 # internal

end

# uinttype
uinttype(::DigitalNet32) = UInt32
uinttype(::DigitalNet64) = UInt64


Base.length(digital_net::DigitalNet32) = digital_net.n
Base.length(digital_net::DigitalNet64) = digital_net.n


const DigitalNet = DigitalNet32

function DigitalNet32(C::Matrix{UInt32}, s::Integer, n::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
#    DigitalNet32{s}(view(z, 1:s), n)
C=reversebits.(C)
t = maximum([ length(string(C[1,i], base=2)) for i in 1:size(C,2) ])
    DigitalNet32{s}(view(C,1:s,:), n,exp2(-t))

end

function DigitalNet64(C::Matrix{UInt64}, s::Integer, n::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
#    DigitalNet32{s}(view(z, 1:s), n)
C=reversebits.(C)
t = maximum([ length(string(C[1,i], base=2)) for i in 1:size(C,2) ])
    DigitalNet64{s}(view(C,1:s,:), n,exp2(-t))

end


function DigitalNet32(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
    DigitalNet32(sobol_Cs,s,2^30)
end

function DigitalNet64_1(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_Cs64,s,2^30)
end

function DigitalNet64(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha2_Bs53,s,2^30)
end

function DigitalNet64_2(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha3_Bs53,s,2^30)
end



#reversebits(n::U) where {U<:Unsigned} = parse(U, reverse(bitstring(n)), base=2)


#@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat}, digital_net::DigitalNet32, k::UInt32)
#    ϕ_k = reversebits(k) * 2.0^(-32) # gray coded radical inverse function in base 2
#    @inbounds for i in 1:length(x)
#        x[i] = ϕ_k * lattice_rule.z[i]
#        x[i] -= floor(x[i]) # mod 1
#    end
#    x
#end


@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat},d::DigitalNet32, k::UInt32)
    k > 0 || return zeros(length(x)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt32,length(x))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:length(x)
                cur[j] ⊻= d.C[j,i]
            end
        end
    end
    for i in 1:length(x)
        x[i] = d.recipid * cur[i]

    end
    return x
end

@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat},d::DigitalNet64, k::UInt64)
    k > 0 || return zeros(length(x)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt64,length(x))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:length(x)
                cur[j] ⊻= d.C[j,i]
            end
        end
    end
    for i in 1:length(x)
        x[i] = d.recipid * cur[i]

    end
    return x
end
