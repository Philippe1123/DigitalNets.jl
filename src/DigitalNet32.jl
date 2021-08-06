


struct DigitalNet32{s} <: AbstractDigitalNets{s}
    C::Matrix{UInt32} # generating matrix
    n::Int64 # max number of points in the Digital net
    recipid::Float64 # internal
    cur::Vector{<:UInt32}
    state::Vector{<:Int32}

end

struct DigitalNet64{s} <: AbstractDigitalNets{s}
    C::Matrix{UInt64} # generating matrix
    n::Int64 # max number of points in the Digital net
    recipid::Float64 # internal
    cur::Vector{<:UInt64}
    state::Vector{<:Int64}

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
cur= zeros(UInt32,s)

    DigitalNet32{s}(view(C,1:s,:), n,exp2(-t),cur,zeros(Int32,1))

end

function DigitalNet64(C::Matrix{UInt64}, s::Integer, n::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
#    DigitalNet32{s}(view(z, 1:s), n)
C=reversebits.(C)
t = maximum([ length(string(C[1,i], base=2)) for i in 1:size(C,2) ])
cur= zeros(UInt64,s)

    DigitalNet64{s}(view(C,1:s,:), n,exp2(-t),cur,zeros(Int64,1))

end


function DigitalNet32(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
    DigitalNet32(sobol_Cs,s,2^20)
end

function DigitalNet64(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_Cs64,s,2^20)
end

function DigitalNet64InterlacedTwo(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha2_Bs53,s,2^20)
end

function DigitalNet64InterlacedThree(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha3_Bs53,s,2^20)
end

function DigitalNet64InterlacedThree64(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha3_Bs64,s,2^20)
end

function DigitalNet64InterlacedFour(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha4_Bs53,s,2^20)
end



function DigitalNet64InterlacedFive(s::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
#    print(@__DIR__())
    DigitalNet64(sobol_alpha5_Bs53,s,2^20)
end#reversebits(n::U) where {U<:Unsigned} = parse(U, reverse(bitstring(n)), base=2)


#@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat}, digital_net::DigitalNet32, k::UInt32)
#    ϕ_k = reversebits(k) * 2.0^(-32) # gray coded radical inverse function in base 2
#    @inbounds for i in 1:length(x)
#        x[i] = ϕ_k * lattice_rule.z[i]
#        x[i] -= floor(x[i]) # mod 1
#    end
#    x
#end


@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat},d::DigitalNet32, k::UInt32)

#t=@elapsed begin

    k > 0 || return zeros(length(x)) # return zero if k == 0
#    @assert k < d.nmax
#println("k ",k)
#println("x ",x)
     k=xor(k,k>>> 1)#grey code reordering
    cur = zeros(UInt32,length(x))
    for i in 1:length(bitstring(k)) #length of row typically 32
        if ( k & (1 << (i - 1) ) ) != 0 #bitwise and check so that only right rows are addressed
#           println("ln ", length(x))
            for j in 1:length(x) #lenght of dim
                cur[j] ⊻= d.C[j,i] # cur = cur ⊻ d.C[j,i]
                x[j] = d.recipid * cur[j]

            #    println(d.C[j,i])
            #    sleep(0.5)
            end
        end
    end

    return x
end

@inline function reset!(d::DigitalNet32)

@inbounds for i in ndims(d)
d.cur[i]=UInt32(0)
end
d.state[1]=0
end

@inline function reset!(d::DigitalNet64)

@inbounds for i in ndims(d)
d.cur[i]=UInt64(0)
end
d.state[1]=0

end


@inline function next!(x::Vector{<:AbstractFloat},d::DigitalNet32, k::UInt32,cur::Vector{<:UInt32})


    k > 0 || return zeros(length(x)),zeros(length(x))
    ctz= trailing_zeros(k)

    @inbounds for j in 1:ndims(d) #lenght of dim
    cur[j] ⊻= d.C[j,ctz+1] # cur = cur ⊻ d.C[j,i]
        x[j] = d.recipid * cur[j]
    end
    return x,cur

end


@inline function getCur(d::DigitalNet64)

    return d.cur
end

@inline function getState(d::DigitalNet64)

    return d.state[1]
end

@inline function getState(d::DigitalNet32)

    return d.state[1]
end

@inline function setCur(d::DigitalNet64,Newcur::Vector{<:UInt64})

for j in 1:ndims(d)
d.cur[j]=Newcur[j]
end
end

@inline function setState(d::DigitalNet64,NewState::Int64)

d.state[1]=NewState
end

@inline function setState(d::DigitalNet32,NewState::Int64)

d.state[1]=NewState
end


@inline function getCur(d::DigitalNet32)

    return d.cur
end

@inline function setCur(d::DigitalNet32,Newcur::Vector{<:UInt32})

for j in 1:ndims(d)
d.cur[j]=Newcur[j]
end
end

@inline function next!(x::Vector{<:AbstractFloat},d::DigitalNet64, k::UInt64,cur::Vector{<:UInt64})

k > 0 || return zeros(length(x)),zeros(length(x))
ctz= trailing_zeros(k)

@inbounds for j in 1:ndims(d) #lenght of dim
    cur[j] ⊻= d.C[j,ctz+1] # cur = cur ⊻ d.C[j,i]
    x[j] = d.recipid * cur[j]
end

return x,cur

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


@inline function unsafe_getcur!(d::DigitalNet64, k::UInt64)
    k > 0 || return zeros(ndims(d)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt64,ndims(d))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:ndims(d)
                cur[j] ⊻= d.C[j,i]
            end
        end
    end
    setCur(d,cur)

end

@inline function unsafe_getcur!(d::DigitalNet32, k::UInt32)
    k > 0 || return zeros(ndims(d)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt32,ndims(d))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:ndims(d)
                cur[j] ⊻= d.C[j,i]
            end
        end
    end
    setCur(d,cur)
end
