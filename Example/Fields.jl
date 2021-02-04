module PlotPoints

using DigitalNets
using PyPlot
using LatticeRules
using Revise
using MultilevelEstimators








function Run()


    correlation_length = 1.0
    smoothness = 3
    covariance_function = CovarianceFunction(2, Matern(correlation_length, smoothness))
    n = 12
   pts = range(0, stop=1, length=n)
    n_kl = 78
    grf = GaussianRandomField(covariance_function, KarhunenLoeve(n_kl), pts, pts)


pt=100

TruncNormal=MultilevelEstimators.TruncatedNormal(0.0,1.0,-1,1)
#TruncNormal=MultilevelEstimators.Uniform(-3,3)

la=LatticeRule(n_kl)
dig_1=dig=DigitalNet64(n_kl)
dig_2=dig=DigitalNet64_1(n_kl)
dig=DigitalNet64_2(n_kl)

dig_1_gauss=zeros(n_kl,1)
for id=1:n_kl

dig_1_gauss[id]=MultilevelEstimators.transform(TruncNormal,dig_1[pt][id])


end
figure()
surf(30e4.+(6e4)*sample(grf,xi=dig_1_gauss))

end
