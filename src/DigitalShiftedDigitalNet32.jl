struct DigitalShiftedDigitalNets32{s, L, V} <: AbstractDigitalNets{s}
    digital_net::L
    Δ::V
end

struct DigitalShiftedDigitalNets64{s, L, V} <: AbstractDigitalNets{s}
    digital_net::L
    Δ::V
end


uinttype(::DigitalShiftedDigitalNets32) = UInt32
uinttype(::DigitalShiftedDigitalNets64) = UInt64

Base.length(d::DigitalShiftedDigitalNets32) = length(d.digital_net)
Base.length(d::DigitalShiftedDigitalNets64) = length(d.digital_net)




DigitalShiftedDigitalNets32(digital_net::DigitalNet32{s}) where s = DigitalShiftedDigitalNets32(digital_net, rand(s)) # specify lattice rule
DigitalShiftedDigitalNets32(s::Integer) = DigitalShiftedDigitalNets32(DigitalNet32(s)) # specify number of dimensions only

DigitalShiftedDigitalNets64(digital_net::DigitalNet64{s}) where s = DigitalShiftedDigitalNets64(digital_net, rand(s)) # specify lattice rule
DigitalShiftedDigitalNets64(s::Integer) = DigitalShiftedDigitalNets64(DigitalNet32(s)) # specify number of dimensions only


function DigitalShiftedDigitalNets32(digital_net::DigitalNet32{s}, Δ::Vector{<:AbstractFloat}) where s
    length(Δ) == ndims(digital_net) || throw(DimensionMismatch("length of the random shift vector must be equal to the number of dimensions of the digital net, expected $(ndims(digital_net)), got $(length(Δ))"))
    all(0 .≤ Δ .≤ 1) || throw(ArgumentError("random shift vector must contain uniformly distributed random numbers"))
    Δ=convert(Array{UInt32},floor.(Δ*typemax(UInt32)))
    DigitalShiftedDigitalNets32{s, typeof(digital_net), typeof(Δ)}(digital_net,Δ )
end

function DigitalShiftedDigitalNets64(digital_net::DigitalNet64{s}, Δ::Vector{<:AbstractFloat}) where s
    length(Δ) == ndims(digital_net) || throw(DimensionMismatch("length of the random shift vector must be equal to the number of dimensions of the digital net, expected $(ndims(digital_net)), got $(length(Δ))"))
    all(0 .≤ Δ .≤ 1) || throw(ArgumentError("random shift vector must contain uniformly distributed random numbers"))
    Δ=convert(Array{UInt64},floor.(Δ*typemax(UInt64)))
    DigitalShiftedDigitalNets64{s, typeof(digital_net), typeof(Δ)}(digital_net,Δ )
end

@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat},d::DigitalShiftedDigitalNets32, k::UInt32)
#    k > 0 || return zeros(length(x)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt32,length(x))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:length(x)
                cur[j] ⊻= d.digital_net.C[j,i]
            end
        end
    end
    for i in 1:length(x)
        x[i] = (cur[i]⊻d.Δ[i])*d.digital_net.recipid
#x[i] = mod((cur[i]+d.Δ[i])*d.digital_net.recipid,1)

    end
    return x
end

@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat},d::DigitalShiftedDigitalNets64, k::UInt64)
#    k > 0 || return zeros(length(x)) # return zero if k == 0
#    @assert k < d.nmax
     k=xor(k,k>>> 1)
    cur = zeros(UInt64,length(x))
    for i in 1:length(bitstring(k))
        if ( k & (1 << (i - 1) ) ) != 0
            for j in 1:length(x)
                cur[j] ⊻= d.digital_net.C[j,i]
            end
        end
    end
    for i in 1:length(x)
        x[i] = (cur[i]⊻d.Δ[i])*d.digital_net.recipid
#        x[i] = mod((cur[i]+d.Δ[i])*d.digital_net.recipid,1)
#       x[i]=cur[i]*d.digital_net.recipid
    end
    return x
end


@inline function next!(x::Vector{<:AbstractFloat},d::DigitalShiftedDigitalNets64, k::UInt64,cur::Vector{<:UInt64})

k > 0 || return zeros(length(x)),zeros(length(x))
ctz= trailing_zeros(k)

@inbounds for j in 1:ndims(d) #lenght of dim
    cur[j] ⊻= d.digital_net.C[j,ctz+1] # cur = cur ⊻ d.C[j,i]



    x[j] = d.digital_net.recipid * cur[j]
end

return x,cur

end


@inline function next!(x::Vector{<:AbstractFloat},d::DigitalShiftedDigitalNets32, k::UInt32,cur::Vector{<:UInt32})

k > 0 || return zeros(length(x)),zeros(length(x))
ctz= trailing_zeros(k)

@inbounds for j in 1:ndims(d) #lenght of dim
    cur[j] ⊻= d.digital_net.C[j,ctz+1] # cur = cur ⊻ d.C[j,i]



    x[j] = d.digital_net.recipid * cur[j]
end

return x,cur

end


@inline function getCur(d::DigitalShiftedDigitalNets64)

    return d.digital_net.cur
end

@inline function getState(d::DigitalShiftedDigitalNets64)

    return d.digital_net.state[1]
end

@inline function getState(d::DigitalShiftedDigitalNets32)

    return d.digital_net.state[1]
end

@inline function setCur(d::DigitalShiftedDigitalNets64,Newcur::Vector{<:UInt64})

for j in 1:ndims(d)
d.digital_net.cur[j]=Newcur[j]
end
end

@inline function setState(d::DigitalShiftedDigitalNets64,NewState::Int64)

d.digital_net.state[1]=NewState
end

@inline function setState(d::DigitalShiftedDigitalNets32,NewState::Int64)

d.digital_net.state[1]=NewState
end


@inline function getCur(d::DigitalShiftedDigitalNets32)

    return d.digital_net.cur
end

@inline function setCur(d::DigitalShiftedDigitalNets32,Newcur::Vector{<:UInt32})

for j in 1:ndims(d)
d.digital_net.cur[j]=Newcur[j]
end
end


@inline function reset!(d::DigitalShiftedDigitalNets32)

@inbounds for i in ndims(d)
d.digital_net.cur[i]=UInt32(0)
end
d.digital_net.state[1]=0
end

@inline function reset!(d::DigitalShiftedDigitalNets64)

@inbounds for i in ndims(d)
d.digital_net.cur[i]=UInt64(0)
end
d.digital_net.state[1]=0

end
