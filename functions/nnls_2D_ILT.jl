using NonNegLeastSquares

"""
	REXSY(spinMC, esp, n1, tm, n2)

	Return steady state magnetizations for REXSY applied to a single spin, over different
    n1, tm, and n2 times. Optional parameters to magnetic field inhomogeniety and gaussian noise.

	Ref: https://pubmed.ncbi.nlm.nih.gov/19894951/

	Inputs:
		spinMC:SpinMC - Multi-compartment spin whose steady state magnetization is being computed
        esp::Real - echo spacing time (in ms) for each n1/n2 loop (in ms)
        n1::Vector{Integer} - Vector of different n1 loop numbers to be simulated (each of duration esp)
		tm::AbstractArray - Vector of different mixing times to be simulated (in ms)
        n2::Vector{Integer} - Vector of different n2 loop numbers to be simulated (each of duration esp)

	Output:
		matrix of complex magnetizations (as an n1 x tm x n2 x noise Array)
"""

function nnls_2D_ILT(REXSY_output::AbstractArray, esp::Int, r1_mesh::AbstractArray,
                    r2_mesh::AbstractArray, t_mesh::AbstractArray)

    tm_length = size(REXSY_output,2)

    t1mat = ones(length(r1_mesh))*t_mesh
    t2mat = ones(length(r2_mesh))*t_mesh

    r1mat = esp*r1_mesh*ones(length(mesh))'
    r2mat = esp*r2_mesh*ones(length(mesh))'

    n1_mat = exp.(-r1mat ./ t1mat)
    n2_mat = exp.(-r2mat ./ t2mat)

    A = kron(n1_mat,n2_mat)

    peak_vec = nonneg_lsq(A,vec(combined_spins[i,:,:]); alg=:fnnls)
    tmp = reshape(peak_vec,(length(t),length(t)))
end
