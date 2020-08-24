

abstract type AbstractDigitalNets{s} <: Random.AbstractRNG end

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





@inline function getpoint(digital_net::AbstractDigitalNets, k::Number) # get the k-th point of the lattice sequence
    0 ≤ k < length(digital_net) || throw(BoundsError(digital_net, k))
    unsafe_getpoint(digital_net, convert(uinttype(digital_net), k))
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



Base.iterate(lattice_rule::AbstractDigitalNets, state=uinttype(lattice_rule)(0)) = state ≥ length(lattice_rule) ? nothing : (getpoint(lattice_rule, state), state + uinttype(lattice_rule)(1))
Base.eltype(::Type{<:AbstractDigitalNets}) = Vector{Float64}

# enable lattice_rule[i] access
Base.getindex(lattice_rule::AbstractDigitalNets, i::Number) = getpoint(lattice_rule, i)
Base.getindex(lattice_rule::AbstractDigitalNets, I) = [lattice_rule[i] for i in I]
Base.firstindex(lattice_rule::AbstractDigitalNets) = uinttype(lattice_rule)(0)
Base.lastindex(lattice_rule::AbstractDigitalNets) = uinttype(lattice_rule)(length(lattice_rule) - 1)
