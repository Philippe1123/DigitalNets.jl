
abstract type DeterministicPoint{s} <: Random.AbstractRNG end

abstract type AbstractDigitalNets{s} <: DeterministicPoint{s} end

# number of dimensions
Base.ndims(::AbstractDigitalNets{s}) where s = s::Int

# size
Base.size(dignet::AbstractDigitalNets) = (length(dignet), )




reversebits(u::UInt32) = begin
    u = ((u >> 1) & 0x55555555) | ((u & 0x55555555) << 1)
    u = ((u >> 2) & 0x33333333) | ((u & 0x33333333) << 2)
    u = ((u >> 4) & 0x0F0F0F0F) | ((u & 0x0F0F0F0F) << 4)
    u = ((u >> 8) & 0x00FF00FF) | ((u & 0x00FF00FF) << 8)
    u = ( u >> 16             ) | ( u               << 16)
end

reversebits(u::UInt64) =begin
u = (UInt64(reversebits(parse(UInt32,string(repr(u))[11:end],base=16) ))<<32) | UInt64(reversebits(parse(UInt32,string(repr(u>>32))[11:end],base=16)));
end


count_trailing_zero_bits(v::UInt32)=begin
       c=32
       v&=-Int32(v)
       if(v!=0) c=c-1 end
       if (v & 0x0000FFFF)!=0 c -= 16 end
       if (v & 0x00FF00FF)!=0 c -= 8 end
       if (v & 0x0F0F0F0F)!=0 c -= 4 end
       if (v & 0x33333333)!=0 c -= 2 end
       if (v & 0x55555555)!=0 c -= 1 end
       return c
       end



@inline function getpoint(digital_net::AbstractDigitalNets, k::Number) # get the k-th point of the lattice sequence
    0 ≤ k < length(digital_net) || throw(BoundsError(digital_net, k))
    unsafe_getpoint(digital_net, convert(uinttype(digital_net), k))
end


@inline function getnextpoint(digital_net::AbstractDigitalNets, k::Number) # get the k-th point of the lattice sequence
    0 ≤ k < length(digital_net) || throw(BoundsError(digital_net, k))
    unsafe_getnextpoint(digital_net, convert(uinttype(digital_net), k))
end

@inline unsafe_getnextpoint(digital_net::AbstractDigitalNets{s}, k::UInt32) where s = begin
    x = Vector{Float64}(undef, s)
    next!(x, digital_net, k) # dispatch to AbstractLatticeRule subtype
end

@inline unsafe_getnextpoint(digital_net::AbstractDigitalNets{s}, k::UInt64) where s = begin
    x = Vector{Float64}(undef, s)
    next!(x, digital_net, k) # dispatch to AbstractLatticeRule subtype
end

# get the k-th point of the sequence without bounds checking for 32 bit integers
@inline unsafe_getpoint(digital_net::AbstractDigitalNets{s}, k::UInt32) where s = begin
    x = Vector{Float64}(undef, s)
    unsafe_getpoint!(x, digital_net, k) # dispatch to AbstractLatticeRule subtype
end

@inline unsafe_getpoint(digital_net::AbstractDigitalNets{s}, k::UInt64) where s = begin
    x = Vector{Float64}(undef, s)
    unsafe_getpoint!(x, digital_net, k) # dispatch to AbstractLatticeRule subtype
end



#Base.iterate(lattice_rule::AbstractDigitalNets, state=uinttype(lattice_rule)(0)) = state ≥ length(lattice_rule) ? nothing : (getpoint(lattice_rule, state), state + uinttype(lattice_rule)(1))
Base.iterate(Dig_net::AbstractDigitalNets, state=uinttype(Dig_net)(0)) = state ≥ length(Dig_net) ? nothing : (getnextpoint(Dig_net, state), state + uinttype(Dig_net)(1))
Base.eltype(::Type{<:AbstractDigitalNets}) = Vector{Float64}

# enable lattice_rule[i] access
Base.getindex(Dig_net::AbstractDigitalNets, i::Number) = getpoint(Dig_net, i)

function Base.getindex(Dig_net::AbstractDigitalNets, I)
      out =[getnextpoint(Dig_net, i) for i in I]
      reset!(Dig_net)
    return   out
end

Base.firstindex(Dig_net::AbstractDigitalNets) = uinttype(Dig_net)(0)
Base.lastindex(Dig_net::AbstractDigitalNets) = uinttype(Dig_net)(length(Dig_net) - 1)
