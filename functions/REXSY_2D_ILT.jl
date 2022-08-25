using NonNegLeastSquares

"""
	peak_mat, Paa, Pbb, Pab, Pba = REXSY_2D_ILT(REXSY_output, esp, r1_mesh, r2_mesh, t_mesh)

	Compute T2-T2 spectral density maps from the output of a multi-compartment
    3-D REXSY scan.

	Inputs:
		REXSY_output:Union{Array{Float64,3},Array{Float64,4}} -
        esp::Real - echo spacing time (in ms) for each n1/n2 loop (in ms)
        r1_mesh::AbstractArray - Mesh encompassing T2 of first compartment (ms)
		r2_mesh::AbstractArray - Mesh encompassing T2 of second compartment (ms)
        t_mesh::AbstractArray - Mesh encompassing T2 of first compartment (ms)

	Output:
		Paa,Pbb,Pba,Pbb - Summed magnitude vectors for exchange between different
        compartments at each tm.

"""

function REXSY_2D_ILT(REXSY_output::Union{Array{Float64,3},Array{Float64,4},Array{Float64,5}},
                    esp::Int, r1_mesh::AbstractArray, r2_mesh::AbstractArray, t_mesh::AbstractArray)

    orig_dim = size(REXSY_output) #Store initial shape

##Create 2D-ILT matrix
    t1mat = ones(length(r1_mesh))*t_mesh'
    t2mat = ones(length(r2_mesh))*t_mesh'

    r1mat = esp*r1_mesh*ones(length(t_mesh))'
    r2mat = esp*r2_mesh*ones(length(t_mesh))'

    n1_mat = exp.(-r1mat ./ t1mat)
    n2_mat = exp.(-r2mat ./ t2mat)

    A = kron(n1_mat,n2_mat)

## Run nnls and reshape outputs to original dimensions
    if ndims(REXSY_output) == 3
        REXSY_output = reshape(REXSY_output, prod(size(REXSY_output)[1:2]), size(REXSY_output,3))
        peak_vec = nonneg_lsq(A,REXSY_output; alg=:fnnls)
        peak_mat = reshape(peak_vec,length(t_mesh),length(t_mesh),orig_dim[3])

    elseif ndims(REXSY_output) == 4
        REXSY_output = reshape(REXSY_output, prod(size(REXSY_output)[1:2]), prod(size(REXSY_output)[3:4]))
        peak_vec = nonneg_lsq(A,REXSY_output; alg=:fnnls)
        peak_mat = reshape(peak_vec,length(t_mesh),length(t_mesh),orig_dim[3],orig_dim[4])

    elseif ndims(REXSY_output) == 5
        REXSY_output = reshape(REXSY_output, prod(size(REXSY_output)[1:2]), prod(size(REXSY_output)[3:5]))
        peak_vec = nonneg_lsq(A,REXSY_output; alg=:fnnls)
        peak_mat = reshape(peak_vec,length(t_mesh),length(t_mesh),orig_dim[3],orig_dim[4],orig_dim[5])

    end


##Find average midpoint between peaks
    tmp = zeros(orig_dim[3])
    for i in 1:orig_dim[3]
        indices = findall(peak_mat[:,:,i,end,end].>0)
        tmp[i] = Int(floor((minimum(indices)[1] + maximum(indices)[1])/2))
    end
    mid = Integer(floor(sum(tmp)/orig_dim[3]))

## Integration of peaks
    Paa = dropdims(sum(peak_mat[1:mid,1:mid,:,:,:];dims=1:2); dims=(1,2))
    Pbb = dropdims(sum(peak_mat[mid+1:length(t_mesh),mid+1:length(t_mesh),:,:,:];dims=1:2); dims=(1,2))
    Pba = dropdims(sum(peak_mat[1:mid,mid+1:length(t_mesh),:,:,:];dims=1:2); dims=(1,2))
    Pab = dropdims(sum(peak_mat[mid+1:length(t_mesh),1:mid,:,:,:];dims=1:2); dims=(1,2))


    return peak_mat,Paa,Pbb,Pab,Pba

end
