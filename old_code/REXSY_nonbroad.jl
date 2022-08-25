using BlochSim
using EasyFit


"""
	REXSY(spin, n1, n2, tm, esp)

	Return steady state magnetization for bssfp applied to a single spin. Use this function if you want to work
	with RF pulses of arbitrary phase cycling factor ϕ.

	Note: An RF phase train with phase cycling factor ϕ is equivalent to a constant-phase RF train with
	a modified off-resonance term.

	Ref: https://pubmed.ncbi.nlm.nih.gov/19894951/

	Inputs:
		spin - Multi-compartment spin whose steady state magnetization is being computed
        esp - echo spacing time (in ms) for each n1/n2 loop (in ms)
        n1 - Vector of different n1 loop numbers to be simulated (each of duration esp)
		n2 - Vector of different n2 loop numbers to be simulated (each of duration esp)
		tm - Vector of different mixing times to be simulated (in ms)

	Output:
		matrix of magnetizations (as a complex number)
"""

tm_length = length(tm)
n1_length = length(n1)
n2_length = length(n2)

spin = SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b),
    (Δf, 0), (1/(kx*Mb), 1/(kx*Ma)))

pspin = deepcopy(spin)
nspin = deepcopy(spin)

pos_spins = fill(spin, n1_length, tm_length, n2_length)
neg_spins = fill(spin, n1_length, tm_length, n2_length)


n1_ctr = 1
tm_ctr = 1
n2_ctr = 1

## Pulse sequence

excite!(pspin, InstantaneousRF(π/2,0)) #give initial π/2 pulse
excite!(nspin, InstantaneousRF(π/2,π)) #give initial π/2 pulse

for i in 1:n1[end] # First CPMG loop n1 times
    freeprecess!(pspin, esp/2) #wait esp/2
    excite!(pspin, InstantaneousRF(π,π/2)) #π pulse
    freeprecess!(pspin, esp/2) #wait esp/2

    freeprecess!(nspin, esp/2) #wait esp/2
    excite!(nspin, InstantaneousRF(π,π/2)) #π pulse
    freeprecess!(nspin, esp/2)

    if i == n1[n1_ctr]
        pos_spins[n1_ctr, :, :] = deepcopy(fill(pspin, 1, tm_length, n2_length))
        neg_spins[n1_ctr, :, :] = deepcopy(fill(nspin, 1, tm_length, n2_length))
        n1_ctr += 1
    end
end

for j in 1:tm_length

    excite!(spin, InstantaneousRF(π/2)) #π/2 pulse
    freeprecess!(spin, tm[tmstep]) #wait tm
    excite!(spin, InstantaneousRF(π/2)) #π/2 pulse

end

for n2step in 1:length(n2)

        for n2counter in 1:n2[n2step] # First CPMG loop n1 times
            freeprecess!(spin, esp/2) #wait esp/2
            excite!(spin, InstantaneousRF(π,π/2)) #π pulse
            freeprecess!(spin, esp/2) #wait esp/2

            if xflip > 0
                pos_urea_spins[tmstep,n1step,n2step] = signal(spin.M[1])
                pos_water_spins[tmstep,n1step,n2step] = signal(spin.M[2])
            else
                neg_urea_spins[tmstep,n1step,n2step] = signal(spin.M[1])
                neg_water_spins[tmstep,n1step,n2step] = signal(spin.M[2])
            end

        end
    end
end


return combined_spins = pos_urea_spins - neg_urea_spins + pos_water_spins - neg_water_spins
