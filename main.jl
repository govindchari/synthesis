using LinearAlgebra

include("synthesis.jl")
include("algorithms.jl")

function f(x)
    return 0.5 * (m*x[1]^2 + l*x[2]^2)
end

function df(x)
    return [m * x[1];l * x[2]]
end

function p(y,r)
    return df(y) - r * (y)
end

n = 6
d = 2

# Smoothness
l = 1000

# Strong convexity
m = 1

# Condition number
k = l / m

x0 = 10 * rand(d,1)
# x0 = [7.659455957006464;0.7941008308411845]
N = 1000

# Run algorithms
A, B, C, rho, P = synthesis(m, l, n, d)
x_syn, iter_syn = synthesized_algorithm(A, B, C, df, N, [x0;zeros(n-d)])
x_nag, iter_nag = nag(df, x0, N, l, m)
x_tm, iter_tm, x1_tm = tm(df, x0, N, l, m)
x_gd, iter_gd = gd(df, x0, N, l, m)

# Compute Residuals
f_res_syn = zeros(N)
f_res_nag = zeros(N)
f_res_tm = zeros(N)
f_res_gd = zeros(N)
x_res_syn = zeros(N)
x_res_nag = zeros(N)
x_res_tm = zeros(N)
x_res_gd = zeros(N)
for i = 1 : N
    global f_res_syn[i] = f(iter_syn[:,i])
    global f_res_nag[i] = f(iter_nag[:,i])
    global f_res_tm[i] = f(iter_tm[:,i])
    global f_res_gd[i] = f(iter_gd[:,i])

    global x_res_syn[i] = norm(iter_syn[:,i])
    global x_res_nag[i] = norm(iter_nag[:,i])
    global x_res_tm[i] = norm(iter_tm[:,i])
    global x_res_gd[i] = norm(iter_gd[:,i])
end

iter_list = LinRange(1,N,10)

rho_tm = 1 - 1 / sqrt(k)
c_tm = (1/rho_tm) * sqrt(dot(x1_tm, x1_tm) - (1/(m * l) * dot(p(x0,m),p(x0,l))))
x1_gd = x0 - (2 / (l+m)) * df(x0)
lmax = maximum(eigvals(P))
lmin = minimum(eigvals(P))
nrmLtil = norm((l-m)*I(n))
nrmC = norm(C)
c_synth = sqrt((lmax/lmin)*(1+nrmLtil^2 * nrmC^2) + (nrmLtil * nrmC^2)/(2 * lmin))

# f_bound_tm = c_tm^2 * (l/2) * rho_tm .^ (2 .* iter_list)
# f_bound_nag = 0.5 * (m + l) * dot(x0,x0) .* exp.(-(iter_list.-1)./sqrt(k))
# f_bound_gd = 0.5 * l * exp.(-(4 .* (iter_list .- 1)) ./ (k + 1)) .* dot(x1_gd,x1_gd)
# f_bound_opt = (m / 2) * ((sqrt(k)-1)/(sqrt(k)+1)).^(2 .* iter_list) .* dot(x0,x0)

x_bound_tm = c_tm * rho_tm .^ iter_list
x_bound_nag = norm(x0) .* sqrt(1 - 1 / sqrt(k)) .^iter_list
x_bound_gd = sqrt.(exp.(-4 .* (iter_list .- 1) ./ (k + 1)) .* dot(x1_gd, x1_gd))
x_bound_opt = sqrt.(((sqrt(k)-1)/(sqrt(k)+1)).^(2 .* iter_list) .* dot(x0,x0))
x_bound_synth = c_synth .* rho .^ iter_list * norm(x0)

rho_tm = (1 - sqrt(m / l))
println("Theoretical Rho TM: ", rho_tm)
println("Theoretical Rho Syn: ", rho)

# figure(dpi=200)
# plot(f_res_syn, label="Synthesized Algorithm", color="blue")
# plot(f_res_nag, label="Nesterov Accelerated Gradient", color="red")
# plot(f_res_tm, label="Triple Momentum", color="green")
# plot(f_res_gd, label="Gradient Descent", color="purple")
# scatter(iter_list, f_bound_gd, color="purple", facecolors="none")
# scatter(iter_list, f_bound_nag, color="red", facecolors="none")
# scatter(iter_list, f_bound_tm, color="green", facecolors="none")
# plot(iter_list, f_bound_opt, color="black", linestyle="--", label="Theoretical Lower Bound")

# grid(true)
# legend()
# title(L"Distance to Optimal for $\kappa=$" * string(k) * " quadratic")
# yscale("log")
# xlabel("Iteration")
# ylabel(L"$f-f^*$")

figure(dpi=200)
plot(x_res_syn, label="Synthesized Algorithm", color="blue")
plot(x_res_nag, label="Nesterov Accelerated Gradient", color="red")
plot(x_res_tm, label="Triple Momentum", color="green")
plot(x_res_gd, label="Gradient Descent", color="purple")
scatter(iter_list, x_bound_gd, color="purple", facecolors="none")
scatter(iter_list, x_bound_nag, color="red", facecolors="none")
scatter(iter_list, x_bound_synth, color="blue", facecolors="none")
scatter(iter_list, x_bound_tm, color="green", facecolors="none")
plot(iter_list, x_bound_opt, color="black", linestyle="--", label="Theoretical Lower Bound")

grid(true)
legend(fontsize="8")
title(L"Distance to Optimal for $\kappa=$" * string(k) * " quadratic")
yscale("log")
xlabel("Iteration")
ylabel(L"$x-x^*$")
ylim([1e-14,1e11])

show()