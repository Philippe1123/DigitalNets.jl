


using GaussianRandomFields
using MultilevelEstimators
using SimpleMultigrid
using Dates
using DigitalNets
using LatticeRules



function Main()
Run(4)
Run(3)
Run(1)
Run(2)
end

function Run(Select::Int64)

    correlation_length = 0.3
    smoothness = 1.0
    cv = CovarianceFunction(2, Matern(correlation_length, smoothness))
    n = 30
   pts = range(0, stop=1, length=n)
    n_kl = 1000
    grf = GaussianRandomField(cv, KarhunenLoeve(n_kl), pts, pts)

    grfs = Vector{typeof(grf)}(undef, 7)
    for i in 1:7
        n = 2^(i+1)
        pts = 1/n:1/n:1-1/n
        grfs[i] = GaussianRandomField(covariance_function, KarhunenLoeve(n_kl), pts, pts)
    end



sample_log = (level, ω) -> sample_lognormal(level, ω, grfs[level + one(level)])

#sample_function = isNL ? isField ? (index,ξ)->Beam_NL_Het_ALL(index, ξ,Lx,Ly


name="Simple_Diffusion_Problem"
timenow = Dates.now()
timenow = Dates.format(timenow, "dd-mm-yyyy-T:HH:MM:SS")
name = string(name,timenow)


nterms=n_kl
#Select=2
if(Select==1)
    pt=DigitalNet64(nterms) #sobol2
    QMCType="Sobol2"
elseif(Select==2)
    pt=DigitalNet64_2(nterms) #Sobol 3
    QMCType="Sobol3"
elseif(Select==3)
    pt=DigitalNet64_1(nterms) #Sobol
    QMCType="Sobol"
elseif(Select==4)
   pt=LatticeRule(nterms) #Lattice
   pt=LatticeRule32("/home/philippe/.julia/packages/LatticeRules/T4oqV/generating_vectors/CKN_250_20.txt",nterms,2^30)
   QMCType="Lattice"
end

name = string(name,QMCType)
println(name)
folder=string("/home/philippe/.julia/dev/DigitalNets/Example/",name)
if(isdir(folder)==false)
mkdir(folder)
else
    rm(folder,recursive=true)
    mkdir(folder)
end

distributions = [TruncatedNormal(0,1,-2,2) for i in 1:n_kl]

#distributions = [MultilevelEstimators.Uniform(0,1) for i in 1:n_kl]

println(distributions)
println(grf)

estimator = Estimator(SL(), QMC(), sample_log, distributions,folder=folder,name=name,point_generator=pt,sample_mul_factor=1.2,nb_of_tols=30,nb_of_shifts=8)
h = run(estimator, 5e-7)
end






function sample_lognormal(level::Level, ω::Vector{<:Real}, grf::GaussianRandomField)


#    ω=ones(length(ω),1)-abs.(2*ω.-1)
#    for id=1:length(ω)
#    ω[id]=transform(TruncatedNormal(0,1,-2,2),ω[id])
#    end

    # solve on finest grid
    z  = sample(grf, xi = ω)
    af = exp.(z)
    Af = elliptic2d(af)
    bf = fill(one(eltype(Af)), size(Af, 1))
    uf = Af\bf
    Qf = uf[length(uf) ÷ 2]
#    println(typeof(Qf))
#GC.gc()

    # compute difference when not on coarsest grid
    dQ = Qf
    if level != Level(0)
        ac = view(af, 2:2:size(af, 1), 2:2:size(af, 2))
        Ac = elliptic2d(ac)
        bc = fill(one(eltype(Af)), size(Ac, 1))
        uc = Ac\bc
        Qc = uc[length(uc) ÷ 2]
        dQ -= Qc
    end
    #GC.gc()



    dQ, Qf
end




Main()
