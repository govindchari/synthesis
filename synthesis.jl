using LinearAlgebra, JuMP, Clarabel, PyPlot

function synthesis(m, l, n, d)
    @assert(n >= 3 * d)

    L = l * I(d)
    M = m * I(d)

    Ltil = L - M
    Ltilpseu = pinv(Ltil)

    r2min = 0.0
    r2max = 1.0

    # Only works for strongly convex functions
    r = 0.0
    Pi = zeros(d, d)

    A = zeros(n,n)
    B = zeros(n,d)
    C = zeros(d,n)
    Pl = zeros(n+d,n+d)
    rho = 100.0

    tol = 1e-10
    while(r2max - r2min > tol)
        r2 = 0.5 * (r2min + r2max)
        lam = r2

        model = Model(Clarabel.Optimizer)
        set_silent(model)
        @variable(model, Ah[1:n,1:n])
        @variable(model, Bh[1:n,1:d])
        @variable(model, Ch[1:d,1:n])
        @variable(model, P[1:(n+d), 1:(n+d)], PSD)

        P11 = P[1:n,1:n]
        P12 = P[1:n,n+1:end]
        P21 = P[n+1:end,1:n]
        P22 = P[n+1:end,n+1:end]

        J1 = zeros(d, n)
        J2 = zeros(d, n)
        J3 = zeros(d, n)
        J1[1:d,1:d] = I(d)
        J2[1:d,d+1:2*d] = I(d)
        J3[1:d,2*d+1:3*d] = I(d)

        lmi = [-r2*P11 -r2*P12 (0.5*J2*Ah - 0.5*lam*Ch)' Ah' (J3*Ah)';
                -r2*P21 (-r2*P22-r*Pi) (0.5*J2*Bh+0.5*lam*Ltilpseu)' Bh' (J3*Bh)';
                0.5*J2*Ah-0.5*lam*Ch (0.5*J2*Bh+0.5*lam*Ltilpseu) (-Ltilpseu-r*Pi) P12' P22';
                Ah Bh P12 -P11 -P12;
                J3*Ah J3*Bh P22 -P21 -P22]

        @constraint(model, Bh .== (Ah - P11) * J1' * inv(M))
        @constraint(model, Ch * J1' .== 1.0 * Array(I(d)))
        @constraint(model, Ch .== J2 * P11)
        @constraint(model, P21 .== J3 * P11)
        @constraint(model, Ch[1:d,1:d] .== 1.0 * Array(I(d)))
        @constraint(model, -lmi >= 0, PSDCone())
        @objective(model, Min, 0)
        optimize!(model)

        if (termination_status(model) == OPTIMAL)
            r2max = r2
            rho = sqrt(r2)
            C = value.(Ch)
            B = value.(P11) \ value.(Bh)
            A = value.(P11) \ value.(Ah) - B * M * C
            Pl = value.(P)  
        else
            r2min = r2
        end
        println("Rho: ", sqrt(r2), " ", termination_status(model))
    end

    return A, B, C, rho, Pl
end
