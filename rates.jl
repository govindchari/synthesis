using PyPlot

N = 50
k = LinRange(1,N,N)

c1 = 1e2
c2 = 1e7

rho1 = 0.8
rho2 = 0.5

figure(dpi=200)
plot(k, c1.*rho1 .^k, label=L"$c=10^2$ $\rho=0.8$", color="green")
plot(k, c1.*rho2 .^k, label=L"$c=10^2$ $\rho=0.5$", color="red")

yscale("log")
legend()
grid(true)
xlabel("Iterations")
ylabel(L"$\|x^*-x_k\|$")
title("Linear Convergence Rates")

figure(dpi=200)
rho1 = 0.8
rho2 = 0.75
N = 100
k2 = LinRange(1,N,N)
plot(k2, c1.*rho1 .^k2, label=L"$c=10^2$ $\rho=0.8$", color="black")
plot(k2, c2.*rho2 .^k2, label=L"$c=10^7$ $\rho=0.75$", color="blue")

yscale("log")
legend()
grid(true)
xlabel("Iterations")
ylabel(L"$\|x^*-x_k\|$")
title("Linear Convergence Rates")

figure(dpi=200)
xx = LinRange(-2,2,500)
plot(xx, 1.0 * xx.^2, label=L"2-smooth ($y=x^2$)", color="blue")
plot(xx, 100 * xx.^2, label=L"200-smooth ($y=100x^2$)", color="red")
plot(xx, abs.(xx), label=L"Not $\lambda$-smooth ($y=|x|$)", color="black")
legend(fontsize="8", loc ="upper left")

ylim((0, 10))
xlabel(L"$x$")
ylabel(L"$y$")
title("Smoothness")
grid(true)

figure(dpi=200)
xx = LinRange(0,4,100)
kk = 10 .^ xx

gd(k) = sqrt(1 - 1 / k)
nag(k) = sqrt(1 - 1 / sqrt(k))
tm(k) = 1 - 1 / sqrt(k)
lb(k) = ((sqrt(k)-1)/(sqrt(k)+1))

plot(kk, gd.(kk), label="GD", color="blue")
plot(kk, nag.(kk), label="NAG", color="red")
plot(kk, tm.(kk), label="TM", color="green")
plot(kk, lb.(kk), label="Theoretical Lower Bound", color="black", linestyle="dashed")

legend(fontsize="8", loc ="lower right")
xscale("log")

xlabel(L"Condition number $\kappa$")
ylabel(L"Convergence rate $\rho$")
title("Convergence rate vs Condition Number")
grid(true)

figure(dpi=200)
k = LinRange(1, 500, 500)
kappa = 1e3

plot(k, c1 * gd(kappa) .^k, label="GD", color="blue")
plot(k, c1 * nag(kappa).^k, label="NAG", color="red")
plot(k, c1 * tm(kappa).^k, label="TM", color="green")
plot(k, c1 * lb(kappa).^k, label="Theoretical Lower Bound", color="black", linestyle="dashed")

title(L"Convergence Upper Bounds for $\kappa=1000$")

legend(fontsize="8", loc ="lower left")
yscale("log")

xlabel("Iteration")
ylabel(L"$\|x^*-x_k\|$")
grid(true)

show()