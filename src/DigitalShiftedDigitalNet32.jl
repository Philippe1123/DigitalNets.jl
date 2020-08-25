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

    end
    return x
end


#function getPoint(r::RandWrapper{s,q,U,M}, k::N where N<:Integer) where {s,q,U<:Unsigned,M}
#    x = getPoint(r.generator, k)
#    return xor.(r.shifts, repeat(convert(Array{UInt32},floor.(typemax(U)*x)),1,size(r.shifts,2))) * r.generator.recipid
#end



#@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat}, Digshift_DigNet32::DigitalShiftedDigitalNets32, k::UInt32)
#    ϕ_k = reversebits(k) * 2.0^(-32) # gray coded radical inverse function in base 2
#    @inbounds for i in 1:length(x)
#        x[i] = ϕ_k * shifted_lattice_rule.lattice_rule.z[i] + shifted_lattice_rule.Δ[i]
#        x[i] -= floor(x[i]) # mod 1
#    end
#    x
#end
