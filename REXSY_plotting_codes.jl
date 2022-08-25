# n1plot = plot(esp*n1times,urea_spins[:,1,1], label="urea, tₘ = 1")
#     plot!(esp*n1times,water_spins[:,1,1], label="water, tₘ = 1")
#     plot!(esp*n1times,urea_spins[:,length(tm),1], label="urea, tₘ = $(tm[end])")
#     plot!(esp*n1times,water_spins[:,length(tm),1], label="water, tₘ = $(tm[end])")
#     title!("REXSY (n₂ = 1)")
#     xlims!(1,esp*n1times[end])
#     xlabel!("n₁ (ms)")
#     ylabel!("Abolsute value of transverse magnetization")
#     yaxis!(:log)
# display(n1plot)
#
#
# n2plot = plot(esp*n2times,urea_spins[1,1,:], label="urea, tₘ = 1")
#     plot!(esp*n2times,water_spins[1,1,:], label="water, tₘ = 1")
#     plot!(esp*n2times,urea_spins[1,length(tm),:], label="urea, tₘ = $(tm[end])")
#     plot!(esp*n2times,water_spins[1,length(tm),:], label="water, tₘ = $(tm[end])")
#     title!("REXSY (n₁ = 1)")
#     xlims!(1,esp*n2times[end])
#     xlabel!("n₂ (ms)")
#     ylabel!("Abolsute value of transverse magnetization")
#     yaxis!(:log)
# display(n2plot)

# tmplot = plot(tm,urea_spins[1,:,1], label="urea, n₁ = n₂ = 1")
#     plot!(tm,water_spins[1,;,1], label="water, n₁ = n₂ = 1")
#     #plot!(tm,urea_spins[length(n1times),:,n2max], label="urea, n₁ = $(n1times[end]), n₂ = $(n2max)")
#     #plot!(tm,water_spins[length(n1times),:,n2max], label="water, n₁ = $(n1times[end]), n₂ = $(n2max)")
#     title!("REXSY (n₂ = 1)")
#     xlims!(1,maximum(tm))
#     xlabel!("tm (ms)")
#     ylabel!("Abolsute value of transverse magnetization")
#     yaxis!(:log)
# display(tmplot)


# peak_plot = plot(tm,Paa, label = "Paa")
#     plot!(tm,Pbb, label = "Pbb")
#     plot!(tm,Pba, label = "Pba")
#     #plot!(tm,Paa+Pab, label = "Paa+Pab")
#     #plot!(tm,Pbb+Pba, label = "Pbb+Pba")
#     #plot!(tm,Pab, label = "Pab")
#     yaxis!(:log)
#     xlabel!("tₘ (ms)")
#     ylabel!("P(tₘ)")
#     title!("Peak Amplitudes")
#     #ylims!(0.001,500)
#     #xlims!(100,tm[end])
#     ylims!(0.01,5)
#     xlims!(100,tm[end])
# display(peak_plot)

## Sanity checks

# test = fitexp(esp*n1times,urea_spins[1,:,1],n=1);
# test.b
#
# test = fitexp(tm,urea_spins[:,1,1],n=1);
# test.b
#
# test = fitexp(esp*n1times,water_spins[1,:,1],n=1);
# test.b
#
# test = fitexp(tm,water_spins[:,1,1],n=1);
# test.b

#SpinMC(1, (Ma, Mb), (1/r1a, 1/r1b), (1/r2a, 1/r2b), (Δf, 0), (1/kab, 1/kba))
