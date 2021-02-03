using DelimitedFiles

using LinearAlgebra

using GaussianRandomFields

using DigitalNets, LatticeRules

using Plots

using Printf

using Statistics

using LaTeXStrings

using SpecialFunctions

function main()

	# Number of grid points
	npt = 128
	h = 1/(npt-1)

	# Number of dimensions
	dim = 100

	# Gaussian random field
	mat = Matern(1, 3)
	cov = CovarianceFunction(1, mat)
	pts = [i/(npt-1) - h/2 for i in 1:npt-1]
	grf = GaussianRandomField(cov, KarhunenLoeve(dim), pts)

	# Plot some realizations of the Gaussian random field
	sps = [exp.(sample(grf, xi=rand(dim) .- 1/2)) for _ in 1:10]
	#sps = [exp.(sample(grf)) for _ in 1:10] # for an actual Gaussian random field 
	plot(pts, sps[1], legend=false, framestyle=:box, xlabel=L"x", ylabel=L"a(x, \omega)")
	for i in 2:length(sps)
		plot!(pts, sps[i])
	end
	savefig("random-fields.pdf")

	# Plot some PDE solutions
	sol = pde_solve.(sps)
	plot(0:h:1, sol[1], legend=false, framestyle=:box, xlabel=L"x", ylabel=L"u(x, \omega)")
	for i in 2:length(sol)
		plot!(0:h:1, sol[i])
	end
	savefig("pde-solutions.pdf")

	# Define function that returns the QoI
	f(x) = begin
		#x = sqrt(2)*erfinv.(2*x .- 1) # transform to normal for actual Gaussian random field
		a = exp.(sample(grf, xi=x .- 1/2))
		u = pde_solve(a)
		m = cld(npt, 2)
		isodd(npt) ? u[m] : (u[m] + u[m+1])/2
	end

	# Log2 of number of integration points
	nip = 16

	# Log2 of number of shifts
	nsh = 3

	# MC
	qpt = [rand(dim) for _ in 1:2^nip] 
	dev = integrate2(f, "Random sampling", nip, 0, qpt)
	plot(2 .^(1:nip), dev, xaxis=:log, yaxis=:log, label="MC", framestyle=:box, legend=:bottomleft,
		 xlabel="number of points", ylabel="root mean square error")

	# QMC, lattice rule
	gen = parse.(UInt32, readlines("exod2_base2_m20_CKN.txt"))
	lat = LatticeRule(gen, dim)
	slt = [ShiftedLatticeRule(lat) for s in 1:2^nsh]
	qpt = [slt[s][n] for n in 0:2^(nip - nsh + 1)-1, s in 1:2^nsh]
	dev = integrate2(f, "Lattice rule", nip - nsh + 1, nsh, qpt)
	plot!(2 .^(nsh:nip), dev, xaxis=:log, yaxis=:log, label="Lattice rule")

	# QMC, digital net
	dgn = DigitalNet32(dim)
	sdn = [DigitalShiftedDigitalNets32(dgn) for s in 1:2^nsh]
	qpt = [sdn[s][n] for n in 0:2^(nip - nsh + 1)-1, s in 1:2^nsh]
	dev = integrate2(f, "Digital net", nip - nsh + 1, nsh, qpt)
	plot!(2 .^(nsh:nip), dev, xaxis=:log, yaxis=:log, label="Digital net")

	# QMC, digital net, order 2
	dgn = DigitalNet64(dim)
	sdn = [DigitalShiftedDigitalNets64(dgn) for s in 1:2^nsh]
	qpt = [sdn[s][n] for n in 0:2^(nip - nsh + 1)-1, s in 1:2^nsh]
	dev = integrate2(f, "Digital net (order 2)", nip - nsh + 1, nsh, qpt)
	plot!(2 .^(nsh:nip), dev, xaxis=:log, yaxis=:log, label="Digital net (order 2)")

	# QMC, digital net, order 3 
	dgn = DigitalNet64_2(dim)
	sdn = [DigitalShiftedDigitalNets64(dgn) for s in 1:2^nsh]
	qpt = [sdn[s][n] for n in 0:2^(nip - nsh + 1)-1, s in 1:2^nsh]
	dev = integrate2(f, "Digital net (order 3)", nip - nsh + 1, nsh, qpt)
	plot!(2 .^(nsh:nip), dev, xaxis=:log, yaxis=:log, label="Digital net (order 3)")
	
	# Plot the result
	savefig("higher-order.pdf")
end

function pde_solve(a)
	n = length(a)
	A = SymTridiagonal(view(a, 1:n-1) + view(a, 2:n), -view(a, 2:n-1))
	x = fill(1/n^2, size(A, 1))
	fct = ldlt(A)
	ldiv!(fct, x)
	prepend!(x, 0)
	append!(x, 0)
	return x
end

function integrate2(fun, name, nip, nsh, qpt)
	dev = zeros(nip)
	qoi = zeros(2^nip, 2^nsh)
	npr = 1
	for ell in 1:nip
		nxt = cld(2^ell, 2) # 1, 1, 2, 4, 8, ...
		for s in 1:2^nsh
			for n in npr:nxt
				qoi[n, s] = fun(qpt[n, s])
			end
		end
		mnr = mean(view(qoi, 1:nxt, :), dims=1)
		mel = mean(mnr)
		sel = nsh == 0 ? std(view(qoi, 1:nxt, :)) / sqrt(nxt) : std(mnr) / sqrt(2^nsh)
		@printf("%s with %d points, mean = %7.5f, std estimator = %13.7e\n", name, nxt, mel, sel)
		dev[ell] = sel
		npr = nxt + 1
	end
	dev
end

main()
