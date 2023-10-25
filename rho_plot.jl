include("synthesis.jl")

N = 50
xx = LinRange(0.01,4,N)
kk = 10 .^ xx

tm(k) = 1 - 1 / sqrt(k)
lb(k) = ((sqrt(k)-1)/(sqrt(k)+1))

tm_rho = tm.(kk)
synth_rho = zeros(N)
for i = 1 : N
    _, _, _, rho, _ = synthesis(1, kk[i], 6, 2)
    synth_rho[i] = rho
end

figure(dpi=200)
plot(kk, tm.(kk), label="Triple Momentum", color="red")
plot(kk, synth_rho, label="Synthesized", color="blue", linestyle="dashed")
plot(kk, lb.(kk), label="Theoretical Lower Bound", color="black", linestyle="dashed")
legend(fontsize="8", loc ="lower right")
xscale("log")

xlabel(L"Condition number $\kappa$")
ylabel(L"Convergence rate $\rho$")
title("Convergence rate vs Condition Number")
grid(true)
show()