function synthesized_algorithm(A, B, C, df, N, x0)
    d, n = size(C)
    x = copy(x0)
    iter = zeros(d, N)
    for i = 1 : N
        x = A * x + B * df(C * x)
        iter[:,i] = C * x
    end
    return C*x, iter
end

function nag(df, x0, N, l, m)
    n = size(x0)[1]
    k = l / m
    b = (sqrt(k) - 1) / (sqrt(k) + 1)
    x = zeros(n, N+1)
    y = zeros(n, N+1)
    x[:,1] = copy(x0)
    y[:,1] = copy(x0)
    for i = 1 : N
        x[:,i+1] = y[:,i] - (1/l) * df(y[:,i])
        y[:,i+1] = (1+b)*x[:,i+1] - b * x[:,i]
    end
    return x[:,end], y
end

function nag_quad(df, x, N, l, m)
    a = (2 / (sqrt(l) + sqrt(m)))^2
    b = ((sqrt(m) - sqrt(l)) / (sqrt(m) + sqrt(l))) ^ 2

    n = size(x)[1]
    iter = zeros(n, N)
    m = zeros(n)
    mprev = zeros(n)
    for i = 1 : N
        m = b * mprev + df(x)
        x = x - a * m
        mprev = m
        iter[:,i] = x
    end
    return x, iter
end

function gd(df, x0, N, l, m)
    a = 2 / (m + l)
    n = size(x0)[1]
    iter = zeros(n, N)
    x = copy(x0)
    for i = 1 : N
        x = x - a * df(x)
        iter[:,i] = x
    end
    return x, iter
end

function tm(df, x, N, l, m)
    n = size(x)[1]
    k = l / m
    rho = 1 - 1 / sqrt(k)
    a = (1 + rho) / l
    b = rho^2 / (2 - rho)
    g = rho^2 / ((1 + rho) * (2 - rho))
    d = rho^2 / (1 - rho^2)

    iter = zeros(n, N + 1)

    xi = zeros(n, N+2)
    xi[:,1] = copy(x)
    xi[:,2] = copy(x)
    y = copy(x)
    iter[:,1] = copy(x)
    xk = zeros(n)

    x1 = zeros(n)

    for i = 2 : N + 1
        xi[:,i+1] = (1 + b) * xi[:,i] - b * xi[:,i-1] - a * df(y)
        y = (1 + g) * xi[:,i+1] - g * xi[:,i]
        xk = (1 + d) * xi[:,i+1] - d * xi[:,i]
        if (i == 2)
            x1 = xk
        end
        iter[:,i] = xk
    end
    return xk, iter, x1
end