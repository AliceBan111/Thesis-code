function VARir(VAR::Dict, VARopt::Dict)
    # =========================================================================
    # Compute impulse responses (IRs) for a VAR model
    # =========================================================================

    nsteps = VARopt[:nsteps]
    impact = VARopt[:impact]
    shut   = VARopt[:shut]
    recurs = VARopt[:recurs]

    Fcomp  = VAR[:Fcomp]
    nvar   = VAR[:nvar]
    nlag   = VAR[:nlag]
    sigma  = VAR[:sigma]

    IR = fill(NaN, nsteps, nvar, nvar)  # (horizon, nvar, nshocks)

    # -------------------- Compute Wold representation ----------------------
    PSI = zeros(nvar, nvar, nsteps)
    VAR[:Fp] = zeros(nvar, nvar, nsteps)

    I_index = VAR[:const] + 1
    for ii in 1:nsteps
        if ii <= nlag
            VAR[:Fp][:, :, ii] = VAR[:F][:, I_index:I_index+nvar-1]
        else
            VAR[:Fp][:, :, ii] = zeros(nvar, nvar)
        end
        I_index += nvar
    end

    PSI[:, :, 1] .= I(nvar)
    for ii in 2:nsteps
        aux = zeros(nvar, nvar)
        for jj in 1:ii-1
            aux += PSI[:, :, ii-jj] * VAR[:Fp][:, :, jj]
        end
        PSI[:, :, ii] .= aux
    end
    VAR[:PSI] = PSI

    # ------------------------ Identification: B matrix ---------------------
    ident = VARopt[:ident]

    B = nothing
    if ident == "short"
        B = cholesky(sigma, :L).L
    elseif ident == "long"
        Finf_big = inv(I(length(Fcomp)) - Fcomp)
        Finf = Finf_big[1:nvar, 1:nvar]
        D = cholesky(Finf * sigma * Finf').L
        B = Finf \ D
    elseif ident == "sign"
        if !haskey(VAR, :B)
            error("B matrix required from SignRestrictions")
        else
            B = VAR[:B]
        end
    elseif ident == "iv"
        up = VAR[:resid][:, 1]
        uq = VAR[:resid][:, 2:end]
        aux, fo, lo = CommonSample(hcat(up, VAR[:IV][VAR[:nlag]+1:end, :]))
        p = aux[:, 1]
        q = uq[end-length(p)+1:end, :]
        Z = aux[:, 2:end]

        FirstStage = OLSmodel(p, Z)
        p_hat = FirstStage[:yhat]

        Biv = zeros(nvar)
        Biv[1, 1] = 1.0
        sqsp = zeros(nvar-1)
        for ii in 2:nvar
            SecondStage = OLSmodel(q[:, ii-1], p_hat)
            Biv[ii, 1] = SecondStage[:beta][2]
            sqsp[ii-1] = SecondStage[:beta][2]
        end

        sigma_b = (1/(size(hcat(p,q),1)-VAR[:ntotcoeff])) * ((hcat(p,q) .- mean(hcat(p,q), dims=1))' * (hcat(p,q) .- mean(hcat(p,q), dims=1)))
        S11 = sigma_b[1,1]
        S21 = sigma_b[2:end, 1]
        S22 = sigma_b[2:end, 2:end]
        Q = sqsp * S11 * sqsp' - (S21 * sqsp' + sqsp * S21') + S22
        sp = sqrt(S11 - (S21 - sqsp * S11)' * (Q \ (S21 - sqsp * S11)))
        Biv .*= sp
        B = zeros(nvar, nvar)
        B[:, 1] .= Biv
        VAR[:FirstStage] = FirstStage
        VAR[:sigma_b] = sigma_b
        VAR[:Biv] = Biv
    elseif ident == "exog"
        B = zeros(nvar, nvar)
        B[:, 1] .= VAR[:F][:, VAR[:const] + VAR[:ncoeff] + 1]
    else
        error("Identification incorrectly specified")
    end

    # -------------------- Compute the impulse responses -------------------
    for mm in 1:nvar
        response = zeros(nvar, nsteps)
        impulse = zeros(nvar)
        if impact == 0
            impulse[mm] = 1.0
        elseif impact == 1
            impulse[mm] = 1.0 / B[mm, mm]
        else
            error("impact must be 0 or 1")
        end
        response[:, 1] .= B * impulse
        if shut != 0
            response[shut, 1] .= 0.0
        end

        if recurs == "wold"
            for kk in 2:nsteps
                response[:, kk] .= PSI[:, :, kk] * B * impulse
            end
        elseif recurs == "comp"
            for kk in 2:nsteps
                FcompN = Fcomp^(kk-1)
                response[:, kk] .= FcompN[1:nvar, 1:nvar] * B * impulse
            end
        end
        IR[:, :, mm] .= response'
    end

    VAR[:B] = B
    return IR, VAR
end
    