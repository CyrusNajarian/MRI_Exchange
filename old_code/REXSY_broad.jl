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

pos_spins = fill(spin, n1_length, tm_length, n2_length)
neg_spins = fill(spin, n1_length, tm_length, n2_length)

## Pulse sequence

 #give initial π/2 pulse
excite.(pos_spins[:,1,1], fill(InstantaneousRF(π/2,0), n1_length, 1, 1)) #give initial π/2 pulse

for xflip in [0,π]
    for n1step in 1:length(n1)
        for n2step in 1:length(n2)
            for tmstep in 1:length(tm)

                for i in 1:n1[n1step] # First CPMG loop n1 times
                    freeprecess!.(pos_spins, esp/2) #wait esp/2
                    excite!(spin, InstantaneousRF(π,π/2)) #π pulse
                    freeprecess!(spin, esp/2) #wait esp/2
                end

                excite!(spin, InstantaneousRF(π/2)) #π/2 pulse
                freeprecess!(spin, tm[tmstep]) #wait tm
                excite!(spin, InstantaneousRF(π/2)) #π/2 pulse

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
    end
end


return combined_spins = pos_urea_spins - neg_urea_spins + pos_water_spins - neg_water_spins
