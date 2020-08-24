struct ShiftedDigitalNets32{s, L, V} <: AbstractDigitalNets{s}
    lattice_rule::L
    Δ::V
end


uinttype(::ShiftedDigitalNets32) = UInt32
uinttype(::ShiftedDigitalNets64) = UInt64





ShiftedDigitalNets32(lattice_rule::LatticeRule32{s}) where s = ShiftedLatticeRule32(lattice_rule, rand(s)) # specify lattice rule
