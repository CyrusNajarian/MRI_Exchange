using NonNegLeastSquares
using Plots

"""
	Paa, Pbb, Pab, Pba, [T2_plots] = REXSY_2D_ILT(REXSY_output, esp, r1_mesh, r2_mesh, t_mesh, [tm,plotting])

	Compute T2-T2 spectral density maps from the output of a multi-compartment
    3-D REXSY scan.

	Inputs:
		REXSY_output:SpinMC - Multi-compartment spin whose steady state magnetization is being computed
        esp::Real - echo spacing time (in ms) for each n1/n2 loop (in ms)
        r1_mesh::AbstractArray - Mesh encompassing T2 of first compartment (ms)
		r2_mesh::AbstractArray - Mesh encompassing T2 of second compartment (ms)
        t_mesh::AbstractArray - Mesh encompassing T2 of first compartment (ms)

    Optional inputs:
        tm::AbstractArray - tm vector used in REXSY scan to label T2-T2 plots
        plotting::Bool - Boolean to output T2-T2 plots (default is 1), requires tm

	Output:
		Paa,Pbb,Pba,Pbb - Summed magnitude vectors for exchange between different
        compartments at each tm.
        T2_plots - Vector containing plots for each tm
"""

function REXSY_2D_ILT(REXSY_output::Array{Float64,3}, esp::Int, r1_mesh::AbstractArray,
                    r2_mesh::AbstractArray, t_mesh::AbstractArray)

    tm_length = size(REXSY_output,3)

    t1mat = ones(length(r1_mesh))*t_mesh'
    t2mat = ones(length(r2_mesh))*t_mesh'

    r1mat = esp*r1_mesh*ones(length(t_mesh))'
    r2mat = esp*r2_mesh*ones(length(t_mesh))'

    n1_mat = exp.(-r1mat ./ t1mat)
    n2_mat = exp.(-r2mat ./ t2mat)

    A = kron(n1_mat,n2_mat)

    Paa, Pab, Pba, Pbb = zeros(tm_length), zeros(tm_length), zeros(tm_length), zeros(tm_length)

    for i = 1:tm_length
        peak_vec = nonneg_lsq(A,vec(REXSY_output[:,:,i]); alg=:fnnls)
        peak_mat = reshape(peak_vec,(length(t_mesh),length(t_mesh)))

        if i == 1
            indices = findall(peak_mat.>0)
            global mid = Int(floor((minimum(indices)[1] + maximum(indices)[1])/2))
        end

        ## Integration of peaks
        Paa[i] = sum(peak_mat[1:mid,1:mid])
        Pbb[i] = sum(peak_mat[mid+1:length(t_mesh),mid+1:length(t_mesh)])
        Pba[i] = sum(peak_mat[1:mid,mid+1:length(t_mesh)])
        Pab[i] = sum(peak_mat[mid+1:length(t_mesh),1:mid])

    end

    return Paa,Pbb,Pab,Pba

    end


function REXSY_2D_ILT(REXSY_output::Array{Float64,3}, esp::Int, r1_mesh::AbstractArray,
                    r2_mesh::AbstractArray, t_mesh::AbstractArray,
                    tm::AbstractArray, plotting::Bool=0)

    tm_length = size(REXSY_output,2)

    t1mat = ones(length(r1_mesh))*t_mesh
    t2mat = ones(length(r2_mesh))*t_mesh

    r1mat = esp*r1_mesh*ones(length(mesh))'
    r2mat = esp*r2_mesh*ones(length(mesh))'

    n1_mat = exp.(-r1mat ./ t1mat)
    n2_mat = exp.(-r2mat ./ t2mat)

    A = kron(n1_mat,n2_mat)

    Paa, Pab, Pba, Pbb = zeros(tm_length), zeros(tm_length), zeros(tm_length), zeros(tm_length)
    T2_plots = Array{Plots.Plot{Plots.GRBackend}}(undef,tm_length)

    for i = 1:tm_length
        peaks = nonneg_lsq(A,vec(combined_spins[i,:,:]); alg=:fnnls)
        tmp = reshape(peaks,(length(t_mesh),length(t_mesh)))

        if i == 1
            indices = findall(tmp.>0)
            global mid = Int(floor((minimum(indices)[1] + maximum(indices)[1])/2))
        end


        T2_plots[i] = plot(heatmap(t_mesh',t_mesh',tmp), type = "heatmap", legend = false);
                            vline!([t[mid]]);
                            hline!([t[mid]]);
                            xlabel!("T̃₂⁽¹⁾ (ms)");
                            ylabel!("T̃₂⁽²⁾ (ms)");
                            title!("T2-T2 spectrum for tₘ = $(tm[i])");

        ## Integration of peaks
        Paa[i] = sum(tmp[1:mid,1:mid])
        Pbb[i] = sum(tmp[mid+1:length(t),mid+1:length(t)])
        Pba[i] = sum(tmp[1:mid,mid+1:length(t)])
        Pab[i] = sum(tmp[mid+1:length(t),1:mid])

    end

return Paa,Pbb,Pab,Pba, T2_plots

end
