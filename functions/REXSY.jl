using BlochSim
using Statistics: mean

"""
	magnetization = REXSY(spinMC, esp, n1, n2, tm, [κ, SNR])

	Return steady state magnetizations for REXSY applied to a single spin, over different
    n1, tm, and n2 times. Optional parameters to magnetic field inhomogeniety and gaussian noise.

	Ref: https://pubmed.ncbi.nlm.nih.gov/19894951/

	Inputs:
		spinMC:SpinMC - Multi-compartment spin whose steady state magnetization is being computed
        esp::Real - echo spacing time (in ms) for each n1/n2 loop (in ms)
        n1::Vector{Integer} - Vector of different n1 loop numbers to be simulated (each of duration esp)
        n2::Vector{Integer} - Vector of different n2 loop numbers to be simulated (each of duration esp)
        tm::AbstractArray - Vector of different mixing times to be simulated (in ms)

    Optional inputs:
        κ:Real - inhomogeniety scaling factor (default is 1)
        SNR::Tuple{Float64,Integer} - Desired SNR of mean output signal and number of generated noisy samples

	Output:
		magnetization = matrix of complex magnetizations (as an n1 x tm x n2 x noise samples)
"""

function REXSY(spinMC::SpinMC, esp::Real, n1::Vector{Int}, n2::Vector{Int}, tm::AbstractArray,
                κ::Real=1)

    n1_length = length(n1)
    tm_length = length(tm)
    n2_length = length(n2)

    magnetization = zeros(ComplexF64, n1_length, n2_length, tm_length)

    for xflip in [0,π]
        for n1step in 1:n1_length
            for tmstep in 1:tm_length
                for n2step in 1:n2_length

                    spin = deepcopy(spinMC)

                    excite!(spin, InstantaneousRF(κ*π/2,xflip)) #give initial π/2 pulse

                    for i in 1:n1[n1step] # First CPMG loop n1 times
                        freeprecess!(spin, esp/2) #wait esp/2
                        excite!(spin, InstantaneousRF(κ*π,π/2)) #π pulse
                        freeprecess!(spin, esp/2) #wait esp/2
                    end

                    excite!(spin, InstantaneousRF(κ*π/2)) #π/2 pulse
                    freeprecess!(spin, tm[tmstep]) #wait tm
                    excite!(spin, InstantaneousRF(κ*π/2)) #π/2 pulse

                    for n2counter in 1:n2[n2step] # First CPMG loop n1 times
                        freeprecess!(spin, esp/2) #wait esp/2
                        excite!(spin, InstantaneousRF(κ*π,π/2)) #π pulse
                        freeprecess!(spin, esp/2) #wait esp/2
                    end

                    if xflip > 0
                        magnetization[n1step,n2step,tmstep] += (signal(spin.M[1]) + signal(spin.M[2]))
                    else
                        magnetization[n1step,n2step,tmstep] -= (signal(spin.M[1]) + signal(spin.M[2]))
                    end

                end
            end
        end
    end
    return magnetization
end


function REXSY(spinMC::SpinMC, esp::Real, n1::Vector{Int}, n2::Vector{Int}, tm::AbstractArray,
                SNR::Tuple{Real,Integer}, κ::Real=1)

    magnetization = REXSY(spinMC, esp, n1, n2, tm, κ)

    n1_length = length(n1)
    n2_length = length(n2)
    tm_length = length(tm)

    noise_length = SNR[2]
    SNR_mag = SNR[1]

    if noise_length == 1
        σ = mean(abs.(magnetization))/(10^(SNR_mag/20))
        magnetization = magnetization + σ * randn(ComplexF64, n1_length, n2_length, tm_length) #add noise
    else
        σ = mean(abs.(magnetization))/(10^(SNR_mag/20))
        magnetization = repeat(magnetization, 1, 1, 1, noise_length)
        magnetization = magnetization + σ * randn(ComplexF64, n1_length, n2_length, tm_length, noise_length) #add noise
    end

    return magnetization
end
