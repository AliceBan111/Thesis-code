using LinearAlgebra
function VARhd(VAR, VARopt)
    # =========================================================================
    # Compute historical decomposition (HD) of the VAR
    # =========================================================================
    sigma   = VAR[:sigma]
    Fcomp   = VAR[:Fcomp]
    constant   = VAR[:const]
    F       = VAR[:Ft]'      
    nvar    = VAR[:nvar]
    nvar_ex = VAR[:nvar_ex]
    nvarXeq = VAR[:nvar] * VAR[:nlag]
    nlag    = VAR[:nlag]
    nlag_ex = VAR[:nlag_ex]
    Y       = VAR[:Y]
    X       = VAR[:X][:, constant+1:nvarXeq+constant]
    nobs    = size(Y, 1)
    
    # Recover B matrix
    B = nothing
    if VARopt[:ident] == "short"
        try
            chol_out = cholesky(sigma)
            B = chol_out.L'
        catch
            error("VCV is not positive definite")
        end
    elseif VARopt[:ident] == "long"
        Finf_big = inv(IMatrix(size(Fcomp,1)) - Fcomp)
        Finf = Finf_big[1:nvar, 1:nvar]
        D = cholesky(Finf * sigma * Finf').L'
        B = Finf \ D
    elseif VARopt[:ident] == "sign"
        if VAR[:B] === nothing
            error("Provide B matrix with SR or SignRestrictions")
        else
            B = VAR[:B]
        end
    elseif VARopt[:ident] == "iv"
        error("Historical decomposition not available for external instruments (iv)")
    else
        error("Identification incorrectly specified")
    end

    # Contribution of each shock
    eps = B \ VAR[:resid]'
    B_big = zeros(nlag*nvar, nvar)
    B_big[1:nvar, :] .= B
    Icomp = hcat(LinearAlgebra.I(nvar), zeros(nvar, (nlag-1)*nvar))
    
    HDshock_big = zeros(nlag*nvar, nobs+1, nvar)
    HDshock = zeros(nvar, nobs+1, nvar)
    
    for j in 1:nvar
        eps_big = zeros(nvar, nobs+1)
        eps_big[j, 2:end] .= eps[j, :]
        for i in 2:nobs+1
            HDshock_big[:, i, j] .= B_big * eps_big[:, i] + Fcomp * HDshock_big[:, i-1, j]
            HDshock[:, i, j] .= Icomp * HDshock_big[:, i, j]
        end
    end

    # Initial value
    HDinit_big = zeros(nlag*nvar, nobs+1)
    HDinit = zeros(nvar, nobs+1)
    HDinit_big[:, 1] .= vec(X[1, :])
    HDinit[:, 1] .= Icomp * HDinit_big[:, 1]
    for i in 2:nobs+1
        HDinit_big[:, i] .= Fcomp * HDinit_big[:, i-1]
        HDinit[:, i] .= Icomp * HDinit_big[:, i]
    end

    # constantant
    HDconstant_big = zeros(nlag*nvar, nobs+1)
    HDconstant = zeros(nvar, nobs+1)
    CC = zeros(nlag*nvar)
    if constant > 0
        CC[1:nvar] .= F[:, 1]
        for i in 2:nobs+1
            HDconstant_big[:, i] .= CC + Fcomp * HDconstant_big[:, i-1]
            HDconstant[:, i] .= Icomp * HDconstant_big[:, i]
        end
    end

    # Linear trend
    HDtrend_big = zeros(nlag*nvar, nobs+1)
    HDtrend = zeros(nvar, nobs+1)
    TT = zeros(nlag*nvar)
    if constant > 1
        TT[1:nvar] .= F[:, 2]
        for i in 2:nobs+1
            HDtrend_big[:, i] .= TT*(i-1) + Fcomp * HDtrend_big[:, i-1]
            HDtrend[:, i] .= Icomp * HDtrend_big[:, i]
        end
    end

    # Quadratic trend
    HDtrend2_big = zeros(nlag*nvar, nobs+1)
    HDtrend2 = zeros(nvar, nobs+1)
    TT2 = zeros(nlag*nvar)
    if constant > 2
        TT2[1:nvar] .= F[:, 3]
        for i in 2:nobs+1
            HDtrend2_big[:, i] .= TT2*(i-1)^2 + Fcomp * HDtrend2_big[:, i-1]
            HDtrend2[:, i] .= Icomp * HDtrend2_big[:, i]
        end
    end

    # Exogenous
    HDexo_big = zeros(nlag*nvar, nobs+1)
    HDexo = zeros(nvar, nobs+1, nvar_ex)
    if nvar_ex > 0
        for ii in 1:nvar_ex
            VARexo = VAR[:X_EX][:, ii]
            EXO = zeros(nlag*nvar)
            EXO[1:nvar] .= F[:, nvar*nlag+constant+ii]
            for i in 2:nobs+1
                HDexo_big[:, i] .= EXO * VARexo[i-1] + Fcomp * HDexo_big[:, i-1]
                HDexo[:, i, ii] .= Icomp * HDexo_big[:, i]
            end
        end
    end

    # All decompositions
    HDendo = HDinit + HDconstant + HDtrend + HDtrend2 + sum(HDexo, dims=3) + sum(HDshock, dims=3)

    # Save and reshape
    HD = Dict()
    HD[:shock] = zeros(nobs+nlag, nvar, nvar)

    println("nobs = ", nobs)
    println("nlag = ", nlag)
    println("HDshock size: ", size(HDshock))
    println("HDshock[i, 2:end, j] length: ", length(HDshock[1, 2:end, 1]))
    println("After vcat length: ", length(vcat(fill(NaN, nlag), vec(HDshock[1, 2:end, 1]))))
    println("HD[:shock] size: ", size(HD[:shock]))

    for i in 1:nvar
        for j in 1:nvar
            HD[:shock][:, j, i] .= vcat(fill(NaN, nlag), vec(HDshock[i, 2:end, j]))
        end
    end
    HD[:init]   = vcat(fill(NaN, nlag-1, nvar), HDinit[:, 1:end]')
    HD[:constant]  = vcat(fill(NaN, nlag, nvar), HDconstant[:, 2:end]')
    HD[:trend]  = vcat(fill(NaN, nlag, nvar), HDtrend[:, 2:end]')
    HD[:trend2] = vcat(fill(NaN, nlag, nvar), HDtrend2[:, 2:end]')
    HD[:exo]   = zeros(nobs+nlag, nvar, nvar_ex)
    for i in 1:nvar_ex
        HD[:exo][:, :, i] .= vcat(fill(NaN, nlag, nvar), HDexo[:, 2:end, i]')
    end
    HD[:endo]   = vcat(fill(NaN, nlag, nvar), HDendo[:, 2:end]')

    VAR[:B] = B
    return HD, VAR
end
