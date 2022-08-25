using PERK: GaussianRFF, krr_train, krr;
using FileIO
using JLD2

#algorithm = "bssfp"
algorithm = "REXSY"
dataset = 99

PERK_urea_results = load("PERK packaged data\\PERK_$(algorithm)_final_data_$dataset.jld2","PERK_urea_results");

H = PERK_urea_results.H
λbest = PERK_urea_results.λbest
kernel = GaussianRFF(H, λbest);

truth = PERK_urea_results.train_data.truth
features = PERK_urea_results.train_data.features
ρbest = PERK_urea_results.ρbest
train = krr_train(truth', features', kernel, ρbest)

spins = PERK_urea_results.test_data.spins
yhat = krr(spins, train, kernel)


r1a = PERK_urea_results.test_data.r1a
r1b = PERK_urea_results.test_data.r1b
kx = PERK_urea_results.test_data.kx

errors = ([r1a, r1b, kx]-yhat) ./ [r1a, r1b, kx]
