using Distributions: Uniform
using CurveFit: exp_fit
using Statistics: std, mean
using BlochSim
using Plots

include("$(pwd())\\functions\\REXSY.jl")
include("$(pwd())\\functions\\REXSY_2D_ILT.jl")

## Define tissue parameter ranges
r1a_range = Uniform(0.0014, 0.0015)
r1b_range = Uniform(0.0012, 0.0013)
r2a_range = Uniform(0.016, 0.018)
r2b_range = Uniform(0.004, 0.006)
Δf_range = Uniform(325, 335)
Ma_range = Uniform(0.175, 0.225)
kx_range = Uniform(0.0009, 0.0011)

## Define tissue phantom properties
Δf = 300
r1a = .001440
r1b = .001240
r2a = .017040
r2b = .005000
Ma = 0.2
Mb = 0.8
kx = .001030
Δω = 0
κ = 1

## Define scan parameters

esp = 2 #echo spacing
tm = range(1, length=20, step = 40) #range of mixing times
t2weightperiods = 10 .^ range(log10(2),log10(250), length=10)
n1 = unique(Integer.(round.(t2weightperiods/esp))) #number of n1 loops
n2 = unique(Integer.(round.(t2weightperiods/esp))) #number of n2 loops

D = length(n1)*length(n2)*length(tm) #length of vectorized samples

Δω_range = Uniform(-5,5) #bulk off-resonance in Hz
κ_range = Uniform(0.9,1.1) #B1+ inhomogeneity scaling factor

test_samples = 5
noise_samples = 3

#Initialize storage arrays
kx_pred = zeros(test_samples,noise_samples)
kx_truth = zeros(test_samples)
REXSY_output = zeros(length(n1),length(n2),length(tm),test_samples,noise_samples)
Paa_ρ1, Pbb_ρ1 = zeros(test_samples,noise_samples), zeros(test_samples,noise_samples)

SNR_range = LinRange(20,150,test_samples)
t_mesh = 1:20:241 #range around two expected T2 values for 2D ILT


for i in 1:test_samples

    #Randomly sample from tissue paramter ranges
    # Δf = rand(Δf_range)
    # r1a = rand(r1a_range)
    # r1b = rand(r1b_range)
    # r2a = rand(r2a_range)
    # r2b = rand(r2b_range)
    # kx = rand(kx_range)
    # Ma = rand(Ma_range)
    # Mb = 1-Ma

    # Δω = rand(Δω_range)
    # κ = rand(κ_range)

    spin = SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf+Δω, Δω), (1/(kx*Mb), 1/(kx*Ma)))
    REXSY_output[:,:,:,i,:] = abs.(REXSY(spin, esp, n1, n2, tm, (SNR_range[i],noise_samples), κ))
    kx_truth[i] = kx
end

peak_mat,Paa,Pbb,_,_ = REXSY_2D_ILT(REXSY_output, esp, n1, n2, t_mesh)

tm_step = 18
plot(heatmap(t_mesh,t_mesh,peak_mat[:,:,tm_step,end,end]), type = "heatmap", legend = false)
    xlabel!("T̃₂⁽¹⁾ (ms)")
    ylabel!("T̃₂⁽²⁾ (ms)")
    title!("T2-T2 spectrum for tₘ = $(tm[tm_step])")

a = 3 #SNR sample
b = 1 #Noise sample
plot(tm, Paa[:,a,b], label = "Paa")
    plot!(tm,Pbb[:,a,b], label = "Pbb")
    yaxis!(:log)
    xlabel!("tₘ (ms)")
    ylabel!("P(tₘ)")
    title!("Peak Amplitudes at SNR = $(SNR_range[a])")
exp_fit(tm,Paa[:,a,b])
exp_fit(tm,Pbb[:,a,b])

for i in 1:test_samples
    for j in 1:noise_samples
        Paa_ρ1[i,j] = exp_fit(tm,Paa[:,i,j])[2]
        Pbb_ρ1[i,j] = exp_fit(tm,Pbb[:,i,j])[2]
    end
end

kx_pred = replace(-Paa_ρ1 - Pbb_ρ1 .- r1a .- r1b, NaN => -0.001)

percent_errors = (kx_pred .- kx_truth)./kx_truth
sample_mean_percent_errors = mean(percent_errors;dims=2)
total_percent_error = mean(sample_mean_percent_errors)

coeff_of_variation = std(kx_pred;dims=2)./mean(kx_pred;dims=2)
total_coeff_of_variation = mean(coeff_of_variation)

RMSE_errors = sqrt.(sum(abs2, kx_pred .- kx_truth, dims=2) / noise_samples)
RMSE_percent = RMSE_errors./kx_truth
total_RMSE_percent = mean(RMSE_percent)



RMSE_plot = plot(RMSE_percent * 100)
    title!("RMSE over scans \n (averaged over 100 added SNR = 50 signals)")
    xlabel!("scan index")
    ylabel!("RMSE (%)")

coeff_of_variation_plot = plot(coeff_of_variation)
    title!("Coefficient of Variation over scans \n (averaged over 100 added SNR = 50 signals)")
    xlabel!("scan index")
    ylabel!("Coefficient of Variation")

percent_error_plot = plot(sample_mean_percent_errors * 100)
    title!("Percent Error over scans \n (averaged over 100 added SNR = 50 signals)")
    xlabel!("scan index")
    ylabel!("Percent Error (%)")
