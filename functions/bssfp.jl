using BlochSim
using LinearAlgebra:I


"""
get_opt_bssfp_design()

Get optimized bssfp scan design.
"""
function get_opt_bssfp_design()

	α = [0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292523868763, 0.17453292519943295,
	0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295,
	0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295,
	0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.17453292519943295, 0.6981317005113884, 0.6981317006955519,
	0.6981317004401415, 0.6981317004087904, 0.6981317004647738, 0.6981317006402664, 0.6981317004105282, 0.6981317006814529,
	0.6981317007977318, 0.698131700672166, 0.6981317005211177, 0.6981317004402081, 0.6981317005337004, 0.6981317007323513,
	 0.6981317006220339, 0.698131700661908, 0.6981317007588633, 0.6981317006197256, 0.6981317007977318]

	ϕ = [-3.0782601051282903, -2.7841473592609107, -2.4805583492526897, -2.170401810560269, -1.8784337163057754, -1.5801616997983006,
	-1.2848868722342726, -0.9796700927081677, -0.6878586068794841, -0.3930651697210255, -0.09201812431693177, 0.20299869308014853,
	0.5049026191594305, 0.7986576653937013, 1.100979632866074, 1.394766988659456, 1.6925010953493698, 1.9887812988842988, 2.2923212550764105,
	 2.591681820212611, 2.8986839403289557, -2.9453449164143937, -2.623345064666762, -2.2712612129898138, -1.9468140362604733,
	 -1.6264681305008342, -1.2946626275565323, -0.9543594974301673, -0.6484225534110698, -0.3143558005704897, 0.02342647082388599,
	  0.32841173707193055, 0.6743548833829934, 1.0102487674197884, 1.3347828270143147, 1.6615784234220567, 1.9769854745053559,
	  2.3268242556200556, 2.6721671168814813, 3.003418542863873]

	return (α, ϕ)
end

"""
	bssfp(spin, α, ϕ, TR, TE)

	Return steady state magnetization for bssfp applied to a single spin. Use this function if you want to work
	with RF pulses of arbitrary phase cycling factor ϕ.

	Note: An RF phase train with phase cycling factor ϕ is equivalent to a constant-phase RF train with
	a modified off-resonance term.

	Ref: W. S. Hinshaw, "Image formation by nuclear magnetic resonance: The sensitive-point method", J. of Appl. Phys., 1976.

	Inputs:
		spin - spin whose steady state magnetization is being computed
		α - Flip angle of RF pulse (in radian)
		ϕ - Phase cycling factor of RF pulse (in radian)
		TR - repetition time (in ms)
		TE - echo time (in ms)

	Output:
		steady-state magnetization (as a complex number)
"""
function bssfp(spin::SpinMC, α::Real, ϕ::Real, TR::Real, TE::Real)


	@assert spin.N == 2 "bSSFP signal (with arbitrary phase ϕ) not implemented yet for N > 2 tissue compartments."

	# create a new spin with modified off-resonance
	off_res_factor = ϕ / (2π * TR * 1e-3)
	new_spin = SpinMC(spin.M0, spin.frac, spin.T1, spin.T2, spin.Δf .- off_res_factor, [1/spin.r[1][2], 1/spin.r[2][1]])

	# Set up matrices
	Rxα, = excite(new_spin, InstantaneousRF(α, 0))
	Afp, bfp = freeprecess(new_spin, TR)
	Ate, bte = freeprecess(new_spin, TE)

	# steady state just before tip down
    R = [Matrix(Rxα.A) zeros(3,3); zeros(3,3) Matrix(Rxα.A)]
    M = (I - Matrix(Afp) * R) \ Vector(bfp)

    # magnetization after tip-down
    M = R * M

    # signal at echo time
    M = Matrix(Ate) * M + Vector(bte)

	return complex(sum(M[1:3:end]), sum(M[2:3:end]))
end
