
const sobol_Cs_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_Cs.col")
const sobol_alpha2_Bs53_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_alpha2_Bs53.col")
const sobol_alpha3_Bs53_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_alpha3_Bs53.col")
const sobol_alpha3_Bs64_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_alpha3_Bs64.col")
const sobol_alpha4_Bs53_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_alpha4_Bs53.col")
const sobol_alpha5_Bs53_file = joinpath(@__DIR__(),"..", "generating_matrices", "sobol_alpha5_Bs53.col")

# re-compile if these files change
include_dependency(sobol_Cs_file)
include_dependency(sobol_alpha2_Bs53_file)
include_dependency(sobol_alpha3_Bs53_file)
include_dependency(sobol_alpha3_Bs64_file)


# default generating vectors
const sobol_Cs64 = readdlm(sobol_Cs_file,UInt64)

const sobol_Cs = readdlm(sobol_Cs_file,UInt32)
const sobol_alpha2_Bs53 = readdlm(sobol_alpha2_Bs53_file,UInt64)
const sobol_alpha3_Bs53 = readdlm(sobol_alpha3_Bs53_file,UInt64)
const sobol_alpha3_Bs64 = readdlm(sobol_alpha3_Bs64_file,UInt64)
const sobol_alpha4_Bs53 = readdlm(sobol_alpha4_Bs53_file,UInt64)
const sobol_alpha5_Bs53 = readdlm(sobol_alpha5_Bs53_file,UInt64)
