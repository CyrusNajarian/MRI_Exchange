
using PERK: GaussianRFF, krr_train, krr;
using Random: randperm;
using FileIO
using JLD2
using Statistics: mean, std

#algorithm = "REXSY"
algorithm = "bssfp"

dataset = 84

train_data = load("Simulated Data\\PERK_$(algorithm)_generated_data_$dataset.jld2", "$(algorithm)_train_data");
test_data = load("Simulated Data\\PERK_$(algorithm)_generated_data_$dataset.jld2", "$(algorithm)_test_data");

features = train_data.features
truth = train_data.truth

Nfit = floor(Int,.8*size(features,1)) # use 70% of the data for fitting, 30% for validation
iperm = randperm(size(features,1))
xfit = features[iperm,:][1:Nfit,:]'
yfit = truth[iperm,:][1:Nfit,:]'
xvalidate = features[iperm,:][(1+Nfit):end,:]'
yvalidate = truth[iperm,:][(1+Nfit):end,:]'

#ρtry = 2. .^ (-29:0.5:-26) #REXSY_3
#λtry = repeat(10 .^ LinRange(0.5, .7, 6),1,size(xfit,1)) #REXSY_3

ρtry = 2. .^ (-220:5:-200)
λtry = repeat(10 .^ LinRange(-1, 1, 10),1,size(xfit,1))
fit = zeros(length(ρtry),size(λtry,1))
H = 1200;

for i in 1:length(ρtry)
    for j in 1:size(λtry,1)
        kernel = GaussianRFF(H,λtry'[:,j]);
        train = krr_train(yfit, xfit, kernel, ρtry[i]) # train with "fit" data
        yhat = krr(xvalidate, train, kernel) # test with "validation" data
        fit[i,j] = sqrt(sum(abs2, yhat[end,:] - yvalidate[end,:]) / sum(abs2, yvalidate[end,:])) #calc NRMSE
    end
    println("$i")
end

best = argmin(fit)
ρbest = ρtry[best[1]];
log2(ρbest)
λbest = λtry[best[2],:];
log10(λbest[1,1])

H = 1200;
kernel = GaussianRFF(H,λbest);
train = krr_train(yfit, xfit, kernel, ρbest) # train with "fit" data

test_samples = size(test_data.features,1)
noise_samples = size(test_data.features,3)

RMSE_errors = zeros(test_samples)
kx_pred = zeros(test_samples ,noise_samples)

for i in 1:test_samples
    yhat = krr(test_data.features[i,:,:], train, kernel) # test with "validation" data
    RMSE_errors[i] = sqrt(sum(abs2, yhat[end,:] .- test_data.truth[i,end])
                / noise_samples) #Calculate RMSE for each test_sample
    kx_pred[i,:] =  yhat[end,:]
end


percent_errors = (kx_pred .- test_data.truth[:,end])./test_data.truth[:,end]
sample_mean_percent_errors = mean(percent_errors;dims=2)
total_percent_error = mean(sample_mean_percent_errors)

coeff_of_variation = std(kx_pred;dims=2)./mean(kx_pred;dims=2)
total_coeff_of_variation = mean(coeff_of_variation)

RMSE_percent = RMSE_errors./test_data.truth[:,end]
total_RMSE_percent = mean(RMSE_percent)

struct PERK_results
    train_data
    test_data
    ρbest
    λbest
    H
    percent_errors
    coeff_of_variation
    RMSE_percent
end

PERK_urea_results = PERK_results(train_data,
                                        test_data,
                                        ρbest,
                                        λbest,
                                        H,
                                        percent_errors,
                                        coeff_of_variation,
                                        RMSE_percent
                                        )

save("PERK packages\\PERK_$(algorithm)_final_data_$dataset.jld2","PERK_urea_results",PERK_urea_results)
