using FileIO: save, load
using JLD2
using Statistics: mean
using Distributions: Uniform,Dirac
using Dates

include("$(pwd())\functions\\bssfp.jl")

## Defien storage structure
struct bssfp_generator_results
    features
    truth
    r1a_range
    r1b_range
    r2a_range
    r2b_range
    Ma_range
    kx_range
    Δf_range
    TR_range
    Δω_range
    κ_range
    α
    ϕ
end

## Urea(a), Water(b) phantom properties
r1a_range = Uniform(0.0014, 0.0020)
r1b_range = Uniform(0.0012, 0.0025)
r2a_range = Uniform(0.015, 0.04)
r2b_range = Uniform(0.00333, 0.02)
Δf_range = Uniform(310, 380)
Ma_range = Uniform(0.05, 0.3)
kx_range = Uniform(0.0008, 0.004)

train_samples = 500000 #multiple of resolution for parameter ranges


## Scan parameters

α, ϕ = get_opt_bssfp_design()
D = length(α)

SNR = 50

# B0 and B1+ information
Δω_range = Uniform(-25,25)  #bulk off-resonance in Hz #uniform -1/2TR to 1/2TR
κ_range = Uniform(0.8,1.2) # B1+ inhomogeneity scaling factor #0.8 to 1.2
TR_range = Dirac(20)

##

bssfp_train_features = Array{Float64}(undef,train_samples,D)
bssfp_train_truth =  Array{Float64}(undef,train_samples,10) #10 is the number of potentially varied parameters

## Method 1
Threads.@threads for i in 1:train_samples

    Δf = rand(Δf_range)
    r1a = rand(r1a_range)
    r1b = rand(r1b_range)
    r2a = rand(r2a_range)
    r2b = rand(r2b_range)
    kx = rand(kx_range)
    Ma = rand(Ma_range)
    Mb = 1-Ma
    Δω = rand(Δω_range)
    κ = rand(κ_range)
    TR = rand(TR_range)
    TE = TR/2

    spin = fill(SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf+Δω, Δω), (1/(kx*Mb), 1/(kx*Ma))), D)
    s = bssfp.(spin, κ * α, ϕ, TR, TE)
    bssfp_train_truth[i,:] = [r1a, r1b, r2a, r2b, Ma, Mb, Δf, κ, Δω, kx]

    σ = mean(abs.(s))/(10^(SNR/20)) #pick σ for SNR of 40
    training_noise = σ * randn(ComplexF64,D)

    bssfp_train_features[i,:] = abs.(s+training_noise)

end

##

bssfp_train_data = bssfp_generator_results(bssfp_train_features,
                                       bssfp_train_truth,
                                       r1a_range,
                                       r1b_range,
                                       r2a_range,
                                       r2b_range,
                                       Ma_range,
                                       kx_range,
                                       Δf_range,
                                       TR_range,
                                       Δω_range,
                                       κ_range,
                                       α,
                                       ϕ)

##


r1a_range = Uniform(0.0014, 0.0015)
r1b_range = Uniform(0.0012, 0.0013)
r2a_range = Uniform(0.016, 0.018)
r2b_range = Uniform(0.004, 0.006)
Δf_range = Uniform(325, 335)
Ma_range = Uniform(0.175, 0.225)
kx_range = Uniform(0.0009, 0.0011)

SNR = 50

Δω_range = Uniform(-25,25)
κ_range = Uniform(0.8,1.2)
TR_range = Dirac(20)

test_samples = 5000
noise_samples = 1000

bssfp_test_features = Array{Float64}(undef,test_samples,D,noise_samples)
bssfp_test_truth =  Array{Float64}(undef,test_samples, 10) #11 is the number of potentially varied parameters

for i in 1:test_samples
    Δf = rand(Δf_range)
    r1a = rand(r1a_range)
    r1b = rand(r1b_range)
    r2a = rand(r2a_range)
    r2b = rand(r2b_range)
    kx = rand(kx_range)
    Ma = rand(Ma_range)
    Mb = 1-Ma
    Δω = rand(Δω_range)
    κ = rand(κ_range)
    TR = rand(TR_range)
    TE = TR/2

    spin = fill(SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf + Δω, Δω), (1/(kx*Mb), 1/(kx*Ma))), D)
    s = bssfp.(spin, κ * α, ϕ, TR, TE)

    σ = sqrt(mean(abs2.(s)))/(10^(SNR/20))
    training_noise = σ * randn(ComplexF64,D,noise_samples)

    bssfp_test_features[i,:,:] = (abs.(repeat(s,1,noise_samples)+training_noise))
    bssfp_test_truth[i,:] = [r1a, r1b, r2a, r2b, Ma, Mb, Δf, κ, Δω, kx]
end

bssfp_test_data = bssfp_generator_results(bssfp_test_features,
                                    bssfp_test_truth,
                                    r1a_range,
                                    r1b_range,
                                    r2a_range,
                                    r2b_range,
                                    Ma_range,
                                    kx_range,
                                    Δf_range,
                                    TR_range,
                                    Δω_range,
                                    κ_range,
                                    α,
                                    ϕ)

save("Simulated Data\\PERK_bssfp_generated_data_$dataset.jld2", "bssfp_train_data", bssfp_train_data,
                "bssfp_test_data", bssfp_test_data)
