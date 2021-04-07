
using DigitalNets
using PyPlot
using LatticeRules
using Revise
using MultilevelEstimators




function Run()

TruncNormal=MultilevelEstimators.TruncatedNormal(0.0,1.0,-2,2)
#TruncNormal=MultilevelEstimators.Uniform(0,1)
#TruncNormal=MultilevelEstimators.Normal()
la=LatticeRule(100)
dig64=DigitalNet64(100)
dig64_1=DigitalNet64_1(100)
dig64_2=DigitalNet64_2(100)

figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(MultilevelEstimators.Normal(),la[id][1]),MultilevelEstimators.transform(MultilevelEstimators.Normal(),la[id][2]),color="blue")
title("Lattice Mapped")


end

figure()
for id=0:2^10-1

scatter(la[id][1],la[id][2],color="blue")
title("Lattice")


end


figure()
for id=0:2^10-1

scatter(dig64_2scatt[id][1],dig64_2[id][2])
title("Sobol_3")


end

figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(TruncNormal,dig64_2[id][1]),MultilevelEstimators.transform(TruncNormal,dig64_2[id][2]))
title("Sobol_3 Gauss")


end


figure()
for id=0:2^10-1

scatter(dig64[id][1],dig64[id][2])
title("Sobol_2")


end


figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(TruncNormal,dig64[id][1]),MultilevelEstimators.transform(TruncNormal,dig64[id][2]))
title("Sobol_2 Gauss")


end

figure()
for id=0:2^10-1

scatter(dig64_1[id][1],dig64_1[id][2])
title("Sobol")
end

figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(TruncNormal,dig64_1[id][1]),MultilevelEstimators.transform(TruncNormal,dig64_1[id][2]))
title("Sobol Gauss")


end


figure()
for id=0:2^10-1
MultilevelEstimators.transform(TruncNormal,la[id][1])
scatter(1-abs(2*MultilevelEstimators.transform(TruncNormal,la[id][1])-1),1-abs(2*MultilevelEstimators.transform(TruncNormal,la[id][2])-1))
title("Lattice Gauss Tent")


end



figure()
for id=0:2^10-1

scatter(1-abs(2*MultilevelEstimators.transform(TruncNormal,dig64_2[id][1])-1),1-abs(2*MultilevelEstimators.transform(TruncNormal,dig64_2[id][2])-1))
title("Sobol 3 Gauss Tent")


end

figure()
for id=0:2^10-1

scatter(1-abs(2*MultilevelEstimators.transform(TruncNormal,dig64_1[id][1])-1),1-abs(2*MultilevelEstimators.transform(TruncNormal,dig64_1[id][2])-1))
title("Sobol Gauss Tent")


end

figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(TruncNormal,1-abs(2*dig64_1[id][1]-1)),MultilevelEstimators.transform(TruncNormal,1-abs(2*dig64_1[id][2]-1)))
title("Sobol Gauss Tent")


end

figure()
for id=0:2^10-1

scatter(MultilevelEstimators.transform(TruncNormal,1-abs(2*dig64_2[id][1]-1)),MultilevelEstimators.transform(TruncNormal,1-abs(2*dig64_2[id][2]-1)))
title("Sobol 3 Gauss Tent")


end

end




Run()
