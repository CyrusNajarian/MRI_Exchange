using Distributions: Uniform, Dirac
using FileIO: save, load
using JLD2
using Dates

include("$(pwd())\\functions\\REXSY.jl")

## Urea(a), Water(b) phantom properties
r1a_range = Uniform(0.0014, 0.0020)
r1b_range = Uniform(0.0012, 0.0025)
r2a_range = Uniform(0.015, 0.04)
r2b_range = Uniform(0.00333, 0.02)
Δf_range = Uniform(310, 380)
Ma_range = Uniform(0.05, 0.4)
kx_range = Uniform(0.0008, 0.004)

#n_samples = 60000 #multiple of resolution for parameter ranges

train_samples = 100000
noise_samples = 100
## Pulse sequence constants
esp = 2 #echo spacing
tm = 10 .^ range(log10(2),log10(250), length=10) #range of mixing times
t2weightperiods = 10 .^ range(log10(2),log10(250), length=10)
n1 = unique(Integer.(round.(t2weightperiods/esp)))
n2 = unique(Integer.(round.(t2weightperiods/esp)))

D = length(n1)*length(n2)*length(tm)
SNR = 50

# B0 and B1+ information
Δω_range = Uniform(-25,25)
κ_range = Uniform(0.8,1.2) # B1+ inhomogeneity scaling factor #0.8 to 1.2

REXSY_train_features = Array{Float64}(undef,D,train_samples,noise_samples)
REXSY_train_truth =  Array{Float64}(undef,train_samples,12)

##Define structure for data storage
struct REXSY_generator_results
    features
    truth
    r1a_range
    r1b_range
    r2a_range
    r2b_range
    Ma_range
    kx_range
    Δf_range
    Δω_range
    κ_range
    esp
    n1times
    tm
    n2times
end
## Pulse sequence

gen_counter = 0

# i,j,k looping through range of initialized r1a/r1b/kx parameters
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

    spin = SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf+Δω, Δω), (1/(kx*Mb), 1/(kx*Ma)))

    REXSY_train_features[:,i,:] = reshape(vec(abs.(REXSY(spin, esp, n1, n2, tm, (SNR, noise_samples), κ))),D,noise_samples)
    REXSY_train_truth[i,:] = [r1a, r1b, r2a, r2b, Ma, Mb, Δf, κ, Δω, kx*Mb, kx*Ma, kx]

    if mod(i,1000) == 0
        global gen_counter += 1000
        println("Completed sample $gen_counter at $(Dates.format(now(), "HH:MM:SS"))")
    end

end

REXSY_train_data = REXSY_generator_results(REXSY_train_features,
                    REXSY_train_truth,
                    r1a_range,
                    r1b_range,
                    r2a_range,
                    r2b_range,
                    Ma_range,
                    kx_range,
                    Δf_range,
                    Δω_range,
                    κ_range,
                    esp,
                    tm,
                    n1,
                    n2)

## Generate urea/water data

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

test_samples = 5000
noise_samples = 100

REXSY_test_features = Array{Float64}(undef,D,test_samples,noise_samples)
REXSY_test_truth =  Array{Float64}(undef,test_samples, 12) #10 is the number of potentially varied parameters

Threads.@threads for i in 1:test_samples

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

    spin = SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf+Δω, Δω), (1/(kx*Mb), 1/(kx*Ma)))

    REXSY_test_features[:,i,:] = reshape(abs.(REXSY(spin, esp, n1, n2, tm, (SNR,noise_samples), κ)),D,noise_samples)
    REXSY_test_truth[i,:] = [r1a, r1b, r2a, r2b, Ma, Mb, Δf, κ, Δω, kx*Mb, kx*Ma, kx]

end

REXSY_test_data = REXSY_generator_results(REXSY_test_features,
                    REXSY_test_truth,
                    r1a_range,
                    r1b_range,
                    r2a_range,
                    r2b_range,
                    Ma_range,
                    kx_range,
                    Δf_range,
                    Δω_range,
                    κ_range,
                    esp,
                    tm,
                    n1,
                    n2)

save("Simulated Data\\PERK_REXSY_generated_data_99.jld2", "REXSY_train_data", REXSY_train_data,
                "REXSY_test_data", REXSY_test_data)
