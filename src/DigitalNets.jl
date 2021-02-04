module DigitalNets
import DelimitedFiles
import Random

using DelimitedFiles
using Random
export AbstractDigitalNets, DigitalNet32,digitalnet_data,DelimitedFiles,reversebits,DigitalNet,getpoint,unsafe_getpoint!, unsafe_getpoint,DigitalNet64,DigitalNet64_2,DigitalShiftedDigitalNets32,DigitalShiftedDigitalNets64
export DigitalNet64_1, DeterministicPoint, uinttype

export Higher_Or_QMC_Example,PlotPoints, RandomDiffusion,Higher_Or_QMC_Example_Normal

export sobol_Cs, sobol_Cs_file,sobol_alpha2_Bs53_file, sobol_alpha2_Bs53,sobol_alpha3_Bs53_file,sobol_alpha3_Bs53, sobol_Cs64

for file in ["common", "DigitalNet32" ,"digitalnet_data","DigitalShiftedDigitalNet32"]
    include(string(file, ".jl"))
end

include("../Example/Higher_Or_QMC_Example.jl")
include("../Example/Higher_Or_QMC_Example_Normal.jl")
include("../Example/PlotPoints.jl")
include("../Example/RandomDiffusion.jl")

end # module
