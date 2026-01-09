using LinearAlgebra
function VARvd(VAR, VARopt)
    # =========================================================================
    # Compute forecast error variance decompositions (VDs) for a VAR model
    # =========================================================================
    nsteps = VARopt[:nsteps]
    ident  = VARopt[:ident]
    Fcomp  = VAR[:Fcomp]
    nlag   = VAR[:nlag]
    nvar   = VAR[:nvar]
    sigma  = VAR[:sigma]
    VD     = zeros(nsteps, nvar, nvar)
    
    # Compute Wold representation
    PSI = zeros(nvar, nvar, nsteps)
    VAR[:Fp] = zeros(nvar, nvar, nsteps)
    I = VAR[:const] + 1
    for kk in 1:nsteps
        if kk <= nlag
            VAR[:Fp][:, :, kk] = VAR[:F][:, I:I+nvar-1]
        else
            VAR[:Fp][:, :, kk] .= 0.0
        end
        I += nvar
    end
    
    # Compute multipliers
    PSI[:, :, 1] .= LinearAlgebra.I(nvar)  # identity matrix
    for kk in 2:nsteps
        aux = zeros(nvar, nvar)
        for jj in 1:kk-1
            aux .+= PSI[:, :, kk-jj] * VAR[:Fp][:, :, jj]
        end
        PSI[:, :, kk] .= aux
    end
    VAR[:PSI] = PSI

    # Identification: Recover B matrix
    B = nothing
    if ident == "short"
        try
            chol_out = cholesky(sigma)
            B = chol_out.L'
        catch
            error("VCV is not positive definite")
        end
    elseif ident == "long"
        Finf_big = inv(IMatrix(size(Fcomp,1)) - Fcomp)
        Finf = Finf_big[1:nvar, 1:nvar]
        D = cholesky(Finf * sigma * Finf').L'
        B = Finf \ D
    elseif ident == "sign"
        if VAR[:B] === nothing
            error("You need to provide the B matrix with SR or SignRestrictions")
        else
            B = VAR[:B]
        end
    elseif ident == "iv"
        error("Forecast error variance decomposition not available with external instruments (iv)")
    else
        error("Identification incorrectly specified: choose 'short', 'long', or 'sign'")
    end

    # Calculate the contribution to the MSE for each shock
    for ii in 1:nvar
        # 1-step ahead variance of forecast error
        MSE = zeros(nvar, nvar, nsteps)
        MSE[:, :, 1] .= sigma
        for nn in 2:nsteps
            MSE[:, :, nn] .= MSE[:, :, nn-1] + PSI[:, :, nn] * sigma * PSI[:, :, nn]'
        end
        
        # Structural forecast error contribution of shock ii
        MSE_shock = zeros(nvar, nvar, nsteps)
        MSE_shock[:, :, 1] .= B[:, ii] * B[:, ii]'
        for nn in 2:nsteps
            MSE_shock[:, :, nn] .= MSE_shock[:, :, nn-1] + PSI[:, :, nn] * MSE_shock[:, :, 1] * PSI[:, :, nn]'
        end
        
        # Compute FEVD
        for nn in 1:nsteps
            for kk in 1:nvar
                VD[nn, ii, kk] = 100 * MSE_shock[kk, kk, nn] / MSE[kk, kk, nn]
            end
        end
    end
    
    VAR[:B] = B
    return VD, VAR
end
