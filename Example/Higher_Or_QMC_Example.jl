module Higher_Or_QMC_Example

using DigitalNets
using LatticeRules

using Statistics
using Random
using DelimitedFiles
using PyPlot


const exod2_base2_m20_CKN_file = joinpath(@__DIR__(),"exod2_base2_m20_CKN.txt")
const exod2_base2_m20_CKN = readdlm(exod2_base2_m20_CKN_file,UInt32)
const shift_file = joinpath(@__DIR__(),"Shifts")
const shifts = readdlm(shift_file,Float64)


f(x,c,b)=reshape(exp.(c *transpose((1:size(x,1)).^(-b)) *reshape(x,size(x,1),size(x,2)*size(x,3))),1,size(x,2),size(x,3))


gexact(s,c,b)= prod(expm1.(c*collect(1:s).^(-b))./(c*(collect(1:s).^(-b))))



function MC(s::Int64,c::Float64,b::Float64)

    println("MC")


N=Int64(1e5)

N_v=2 .^ collect(4:20)
N_v=Int64.(ceil.(2 .^ collect(4:20)/3))
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
t= @elapsed begin
end
t=0
for p=1:length(N_v)
t= @elapsed begin
N=N_v[p]

G=f(rand(s,N),c,b)

Mc_Q=mean(G)
Mc_std=std(G)/sqrt(N)
Matrx_std[p]=Mc_std
exact=gexact(s,c,b)
err=abs(Mc_Q-exact)
end
println("MC_Q ", Mc_Q," error= ",err," std= ", Mc_std," N= ", N)
tvec[p]=t

#figure(1)
#loglog(N_v,Matrx_std)
end
println(tvec)
println(Matrx_std)
return Matrx_std,N_v,tvec

end


function QMC_Lattice(s::Int64,c::Float64,b::Float64)

println("QMC Lattice Rule No Shifting")

N=Int64(2^22)

Lattice=LatticeRule(vec(exod2_base2_m20_CKN),s)
Matrix=zeros(s,N+2)
for id=0:length(Lattice[0:N])
Matrix[:,id+1]=Lattice[id]

end

G=f(Matrix,c,b)

Mc_Q=mean(G)
Mc_std=std(G)/sqrt(N)
exact=gexact(s,c,b)
err=abs(Mc_Q-exact)
println("MC_Q ", Mc_Q," error= ",err," std= ", Mc_std," N= ", N)


end

function QMC_Lattice_Shift(s::Int64,c::Float64,b::Float64)

println("QMC Lattice Rule WITH Shifting")


N=Int64(2^22)

N_v=2 .^ collect(4:16)*8
#N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
for p=1:length(N_v)
    t =@elapsed begin
    M=2^3

N=Int64(N_v[p]/M)



Lattice=LatticeRule(vec(exod2_base2_m20_CKN),s)
Matrix=zeros(s,N,M)
#println(Matrix)
#println(Lattice)
shiftLat=ShiftedLatticeRule(Lattice)
#println(shifts)
#println(size(shifts))
for Nshift=1:M
##println(shifts[:,Nshift])
#shiftLat=ShiftedLatticeRule(Lattice,shifts[:,Nshift])
shiftLat=ShiftedLatticeRule(Lattice)
for id=1:length(shiftLat[0:N])-1
Matrix[:,id,Nshift]=shiftLat[id-1]

end
end

G=f(Matrix,c,b)


QMc_R=mean(G,dims=2)
QMc_Q=mean(QMc_R,dims=3)

QMc_std=std(QMc_R)/sqrt(M)
Matrx_std[p]=QMc_std


exact=gexact(s,c,b)
err=abs(QMc_Q[1]-exact)
println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)

end
tvec[p]=t
end
return Matrx_std,N_v,tvec

end




# function QMC_Lattice_Shift_with_TENT_First(s::Int64,c::Float64,b::Float64)
#
#     println("QMC Lattice Rule WITH TENT (first tent then shift)")
#
#
# N=Int64(2^22)
#
# N_v=2 .^ collect(4:16)*8
# #N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
# Matrx_std=zeros(length(N_v),1)
# tvec=zeros(length(N_v),1)
# for p=1:length(N_v)
#     t =@elapsed begin
#     M=2^3
#
# N=Int64(N_v[p]/M)
#
#
#
# Lattice=LatticeRuleTent(vec(exod2_base2_m20_CKN),s)
# Matrix=zeros(s,N,M)
# #println(Matrix)
# #println(Lattice)
# shiftLat=ShiftedLatticeRuleTent(Lattice)
# #println(shifts)
# #println(size(shifts))
# for Nshift=1:M
# ##println(shifts[:,Nshift])
# #shiftLat=ShiftedLatticeRule(Lattice,shifts[:,Nshift])
# shiftLat=ShiftedLatticeRuleTent(Lattice)
# for id=1:length(shiftLat[0:N])-1
# Matrix[:,id,Nshift]=shiftLat[id-1]
#
# end
# end
#
# G=f(Matrix,c,b)
#
#
# QMc_R=mean(G,dims=2)
# QMc_Q=mean(QMc_R,dims=3)
#
# QMc_std=std(QMc_R)/sqrt(M)
# Matrx_std[p]=QMc_std
#
#
# exact=gexact(s,c,b)
# err=abs(QMc_Q[1]-exact)
# println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)
#
# end
# tvec[p]=t
# end
# return Matrx_std,N_v,tvec
#
# end

function QMC_Lattice_Shift_with_TENT_Last(s::Int64,c::Float64,b::Float64)

println("QMC Lattice Rule WITH TENT (first shift then tent)")


N=Int64(2^22)

N_v=2 .^ collect(4:16)*8
#N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
for p=1:length(N_v)
    t =@elapsed begin
    M=2^3

N=Int64(N_v[p]/M)



Lattice=LatticeRule(vec(exod2_base2_m20_CKN),s)
Matrix=zeros(s,N,M)
#println(Matrix)
#println(Lattice)
shiftLat=ShiftedLatticeRule(Lattice)
#println(shifts)
#println(size(shifts))
for Nshift=1:M
##println(shifts[:,Nshift])
#shiftLat=ShiftedLatticeRule(Lattice,shifts[:,Nshift])
shiftLat=ShiftedLatticeRule(Lattice)
for id=1:length(shiftLat[0:N])-1
x=shiftLat[id-1]
out=ones(length(x),1)-abs.(2*x.-1)

Matrix[:,id,Nshift]=out

end
end

G=f(Matrix,c,b)


QMc_R=mean(G,dims=2)
QMc_Q=mean(QMc_R,dims=3)

QMc_std=std(QMc_R)/sqrt(M)
Matrx_std[p]=QMc_std


exact=gexact(s,c,b)
err=abs(QMc_Q[1]-exact)
println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)

end
tvec[p]=t
end
return Matrx_std,N_v,tvec

end

function QMC_Digital_DShift(s::Int64,c::Float64,b::Float64)

println("QMC Digital Net WITH DShifting")

N=Int64(2^22)

N_v=2 .^ collect(4:16)*8
#N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
t= @elapsed begin
end
t=0
for p=1:length(N_v)
    t=0
    t =@elapsed begin
M=2^3
N=Int64(N_v[p]/M)

DigitalNet=DigitalNet32(s)
Matrix=zeros(s,N,M)
#println(Matrix)
for Nshift=1:M
#println(shifts[:,Nshift])
shiftNet=DigitalShiftedDigitalNets32(DigitalNet)
#shiftLat=ShiftedLatticeRule(Lattice)

for id=1:length(shiftNet[0:N])-1
Matrix[:,id,Nshift]=shiftNet[id-1]

end
end

G=f(Matrix,c,b)
#println(Matrix)
#println(G)
#println(size(G,1))
#println(size(G,2))
#println(size(G,3))

QMc_R=mean(G,dims=2)
QMc_Q=mean(QMc_R,dims=3)

QMc_std=std(QMc_R)/sqrt(M)
Matrx_std[p]=QMc_std

#println(QMc_std[1])
#println(QMc_Q[1])

exact=gexact(s,c,b)
err=abs(QMc_Q[1]-exact)
println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)
end
tvec[p]=t
end
return Matrx_std,N_v,tvec


end



function QMC_Digital_DShift64(s::Int64,c::Float64,b::Float64)

println("QMC Digital Net WITH DShifting")


N=Int64(2^22)

N_v=2 .^ collect(4:16)*8
#N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
t= @elapsed begin
end
t=0
for p=1:length(N_v)
    t=0
    t =@elapsed begin
    M=2^3
N=Int64(N_v[p]/M)

DigitalNet=DigitalNet64InterlacedTwo(s)
Matrix=zeros(s,N,M)
#println(Matrix)
for Nshift=1:M
#println(shifts[:,Nshift])
shiftNet=DigitalShiftedDigitalNets64(DigitalNet)
#shiftLat=ShiftedLatticeRule(Lattice)

for id=1:length(shiftNet[0:N])-1
Matrix[:,id,Nshift]=shiftNet[id-1]

end
end

G=f(Matrix,c,b)
#println(Matrix)
#println(G)
#println(size(G,1))
#println(size(G,2))
#println(size(G,3))

QMc_R=mean(G,dims=2)
QMc_Q=mean(QMc_R,dims=3)

QMc_std=std(QMc_R)/sqrt(M)
Matrx_std[p]=QMc_std

#println(QMc_std[1])
#println(QMc_Q[1])

exact=gexact(s,c,b)
err=abs(QMc_Q[1]-exact)
println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)

end
tvec[p]=t
end
return Matrx_std,N_v,tvec

end


function QMC_Digital_DShift64_2(s::Int64,c::Float64,b::Float64)

println("QMC Digital Net WITH DShifting")


N=Int64(2^22)

N_v=2 .^ collect(4:16)*8
#N_v=Int64.(ceil.(2 .^ collect(4:20)/3))*8
Matrx_std=zeros(length(N_v),1)
tvec=zeros(length(N_v),1)
t= @elapsed begin
end
t=0
for p=1:length(N_v)
    t=0
    t =@elapsed begin
    M=2^3
N=Int64(N_v[p]/M)

DigitalNet=DigitalNet64InterlacedThree(s)
Matrix=zeros(s,N,M)
#println(Matrix)
for Nshift=1:M
#println(shifts[:,Nshift])
shiftNet=DigitalShiftedDigitalNets64(DigitalNet)
#shiftLat=ShiftedLatticeRule(Lattice)

for id=1:length(shiftNet[0:N])-1
Matrix[:,id,Nshift]=shiftNet[id-1]

end
end

G=f(Matrix,c,b)
#println(Matrix)
#println(G)
#println(size(G,1))
#println(size(G,2))
#println(size(G,3))

QMc_R=mean(G,dims=2)
QMc_Q=mean(QMc_R,dims=3)

QMc_std=std(QMc_R)/sqrt(M)
Matrx_std[p]=QMc_std

#println(QMc_std[1])
#println(QMc_Q[1])

exact=gexact(s,c,b)
err=abs(QMc_Q[1]-exact)
println("QMc_Q ", QMc_Q," error= ",err," std= ", QMc_std," N= ", N_v[p], "Shifts= ",M , "Samples= ", N)

end
tvec[p]=t
end
return Matrx_std,N_v,tvec

end


function Run()

s=Int64(100)
b=4.
c=2.
MC_std,N_MC,t1 = MC(s,c,b)
QMC_lat_std,N_QMC_lat,t2 =QMC_Lattice_Shift(s,c,b)
QMC_DigDS_std,N_QMC_DigDS,t3 =QMC_Digital_DShift(s,c,b)
QMC_DigDS64_std,N_QMC_DigDS64,t4 =QMC_Digital_DShift64(s,c,b)
QMC_DigDS64_std_2,N_QMC_DigDS64_2,t5 =QMC_Digital_DShift64_2(s,c,b)
# QMC_lat_TENT_std,N_QMC_lat_TENT,t6=QMC_Lattice_Shift_with_TENT_First(s,c,b)
QMC_lat_TENT_std_Last,N_QMC_lat_TENT_Last,t6=QMC_Lattice_Shift_with_TENT_Last(s,c,b)
println(t1)
println(t2)
println(t3)
println(t4)
println(t5)
println(QMC_DigDS64_std_2)

figure()
vec=[10; 10^5]
str=string("b = ",b," ","c = ",c)
loglog(N_MC,MC_std,"-o",label="MC")
loglog(N_QMC_lat,QMC_lat_std,"-o",label="Lattice")
loglog(N_QMC_DigDS,QMC_DigDS_std,"-o",label="Net")
loglog(N_QMC_DigDS64,QMC_DigDS64_std,"-o",label="Net64")
loglog(N_QMC_DigDS64_2,QMC_DigDS64_std_2,"-o",label="Net64")
# loglog(N_QMC_lat_TENT,QMC_lat_TENT_std,"-o",label="Lattice with tent First")
loglog(N_QMC_lat_TENT_Last,QMC_lat_TENT_std_Last,"-o",label="Lattice with tent Last")

loglog(vec,vec.^(-1/2),label=1/2)
loglog(vec,vec.^(-1),label=1/2)
loglog(vec,vec.^(-1.5),label=1/2)
loglog(vec,vec.^(-2),label=1/2)
title(str)

legend(("MC","Lattice (exo2_base2_m20)","sobol_Cs","sobol_alpha2_BS","sobol_alpha3_BS","Lattice with tent Last","O(n^-1/2)","O(n^-1)","O(n^-1.5)","O(n^-2)"))
xlabel("Number of samples n")
ylabel("root mean square error")

figure()
loglog(MC_std,t1,"-o",label="MC")
loglog(QMC_lat_std,t2,"-o",label="Lattice")
loglog(QMC_DigDS_std,t3,"-o",label="Net")
loglog(QMC_DigDS64_std,t4,"-o",label="Net64")
loglog(QMC_DigDS64_std_2,t5,"-o",label="Net64")

legend(("MC","Lattice (exo2_base2_m20)","sobol_Cs","sobol_alpha2_BS","sobol_alpha3_BS"))


end

end
