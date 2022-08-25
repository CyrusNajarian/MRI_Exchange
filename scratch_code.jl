using Plots
using MIRTjim
using NIfTI

ni = zeros(128,128,40)


for i in 1:21
  ni[:,:,i] = niread("1130_tssfp_altered01.nii\\vol_image0$(string(i,pad=2))echo001.nii")
end

for i in 1:19
  ni[:,:,i+21] = niread("1128_tssfp_altered01.nii\\vol_image0$(string(i,pad=2))echo001.nii")
end
rectangle(a) = Shape([(a+1,a+1),(a+1,a-1),(a-1,a-1),(a-1,a+1)])

x1 = 50
x2 = 55
x3 = 60
x4 = 65
x5 = 70
x6 = 75
x7 = 80

scatter(ni[x1,x1,:], title = "pixel ($x1,$x1)")
scatter(ni[x2,x2,:], title = "pixel ($x2,$x2)")
scatter(ni[x3,x3,:], title = "pixel ($x3,$x3)")
scatter(ni[x4,x4,:], title = "pixel ($x4,$x4)")
scatter(ni[x5,x5,:], title = "pixel ($x5,$x5)")
scatter(ni[x6,x6,:], title = "pixel ($x6,$x6)")
scatter(ni[x7,x7,:], title = "pixel ($x7,$x7)")

jim(ni[:,:,1], title = "scan 1")
      plot!(rectangle(x1), opacity=.5)
      plot!(rectangle(x2), opacity=.5)
      plot!(rectangle(x3), opacity=.5)
      plot!(rectangle(x4), opacity=.5)
      plot!(rectangle(x5), opacity=.5)
      plot!(rectangle(x6), opacity=.5)
      plot!(rectangle(x7), opacity=.5)


plot(jim(ni[:,:,1], title = "scan 1"),
      jim(ni[:,:,4], title = "scan 4"),
      jim(ni[:,:,5], title = "scan 5"),
      jim(ni[:,:,6], title = "scan 6"),
      jim(ni[:,:,7], title = "scan 7"),
      jim(ni[:,:,10], title = "scan 10"),
      jim(ni[:,:,13], title = "scan 13"),
      jim(ni[:,:,16], title = "scan 16"),
      jim(ni[:,:,19], title = "scan 19"),
      showaxis = false, xticks = [], yticks = [], colorbar = false)

# include("$(pwd())\\Pulse sequence code\\bssfp.jl")
#
#
# Δf = 330
# r1a = .001440
# r1b = .001240
# r2a = .017040
# r2b = .005000
# Ma = 0.2
# Mb = 0.8
# kx = .001030
#
# TR_samples = 10000
# noise_samples = 20000
#
# α, ϕ = get_opt_bssfp_design()
# D = length(α)
# TR_range = LinRange(2,25,TR_samples)
# SNR = 50
# Δω = 0
# κ = 1
#
# bssfp_TR_features = Array{Float64}(undef,TR_samples,D,noise_samples)
# bssfp_TR_truth =  Array{Float64}(undef,TR_samples,10)
#
# Threads.@threads for i in 1:TR_samples
#
#     TR = TR_range[i]
#     TE = TR/2
#
#     spin = fill(SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf + Δω, Δω), (1/(kx*Mb), 1/(kx*Ma))), D)
#     s = bssfp.(spin, κ * α, ϕ, TR, TE)
#
#     σ = sqrt(mean(abs2.(s)))/(10^(SNR/20))
#     training_noise = σ * randn(ComplexF64,D,noise_samples)
#
#     bssfp_TR_features[i,:,:] = (abs.(repeat(s,1,noise_samples)+training_noise))
#     bssfp_TR_truth[i,:] = [r1a, r1b, r2a, r2b, Ma, Mb, Δf, κ, Δω, kx]
#
# end
#
# RMSE_errors = zeros(TR_samples)
# kx_pred = zeros(TR_samples ,noise_samples)
#
# for i in 1:TR_samples
#     yhat = krr(bssfp_TR_features[i,:,:], train, kernel) # test with "validation" data
#     RMSE_errors[i] = sqrt(sum(abs2, yhat[end,:] .- bssfp_TR_truth[i,end])
#                 / noise_samples) #Calculate RMSE for each test_sample
#     kx_pred[i,:] =  yhat[end,:]
# end
#
# percent_errors = (kx_pred .- bssfp_TR_truth[:,end])./bssfp_TR_truth[:,end]
# sample_mean_percent_errors = mean(percent_errors;dims=2)
# total_percent_error = mean(sample_mean_percent_errors)
#
# coeff_of_variation = std(kx_pred;dims=2)./mean(kx_pred;dims=2)
# total_coeff_of_variation = mean(coeff_of_variation)
#
# RMSE_percent = RMSE_errors./bssfp_TR_truth[:,end]
# total_RMSE_percent = mean(RMSE_percent)
#
# using Plots
#
# RMSE_plot = plot(TR_range,RMSE_percent * 100)
#     title!("TR vs RMSE \n (averaged over 4000 added SNR = 50 signals)")
#     xlabel!("TR (ms)")
#     ylabel!("RMSE (%)")
#
# coeff_of_variation_plot = plot(TR_range,coeff_of_variation)
#     title!("TR vs Coefficient of Variation \n (averaged over 4000 added SNR = 50 signals)")
#     xlabel!("TR (ms)")
#     ylabel!("Coefficient of Variation")
#
# percent_error_plot = plot(TR_range,sample_mean_percent_errors * 100)
#     title!("TR vs Percent Error \n (averaged over 4000 added SNR = 50 signals)")
#     xlabel!("TR (ms)")
#     ylabel!("Percent Error (%)")
