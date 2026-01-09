function VARmodel(ENDO::Matrix{Float64}, nlag::Int, constant::Int=1, EXOG::Union{Matrix{Float64}, Nothing}=nothing, nlag_ex::Int=0)
    # =======================================================================
    # Perform VAR with OLS 
    # =======================================================================
    
    nobs, nvar = size(ENDO)

    VARopt = VARoption()
    
    VAR = Dict()
    VAR[:ENDO] = ENDO
    VAR[:nlag] = nlag
    VAR[:const] = constant
    
    # Check if there are exogenous variables 
    if EXOG !== nothing
        nobs2, nvar_ex = size(EXOG)
        # Check that ENDO and EXOG are conformable
        if nobs2 != nobs
            error("var: nobs in EXOG-matrix not the same as y-matrix")
        end
        VAR[:EXOG] = EXOG
    else
        nvar_ex = 0
        nlag_ex = 0
        VAR[:EXOG] = nothing
    end
    
    # Save some parameters and create data matrices
    nobse = nobs - max(nlag, nlag_ex)
    VAR[:nobs] = nobse
    VAR[:nvar] = nvar
    VAR[:nvar_ex] = nvar_ex    
    VAR[:nlag] = nlag
    VAR[:nlag_ex] = nlag_ex
    ncoeff = nvar * nlag
    VAR[:ncoeff] = ncoeff
    ncoeff_ex = nvar_ex * (nlag_ex + 1)
    ntotcoeff = ncoeff + ncoeff_ex + constant
    VAR[:ntotcoeff] = ntotcoeff
    VAR[:const] = constant
    
    # Create independent vector and lagged dependent matrix
    Y, X = VARmakexy(ENDO, nlag, constant)
    
    # Create (lagged) exogenous matrix
    if nvar_ex > 0
        X_EX = VARmakelags(EXOG, nlag_ex)
        if nlag == nlag_ex
            X = [X X_EX]
        elseif nlag > nlag_ex
            diff = nlag - nlag_ex
            X_EX = X_EX[diff+1:end, :]
            X = [X X_EX]
        elseif nlag < nlag_ex
            diff = nlag_ex - nlag
            Y = Y[diff+1:end, :]
            X = [X[diff+1:end, :] X_EX]
        end
    end
    
    # OLS estimation equation by equation
    for j in 1:nvar
        Yvec = Y[:, j]
        OLSout = OLSmodel(Yvec, X, 0)  
        
        VAR[Symbol("eq$j")] = Dict(
            :beta => OLSout[:beta],
            :tstat => OLSout[:tstat],  
            :bstd => OLSout[:bstd],
            :tprob => OLSout[:tprob],
            :resid => OLSout[:resid],
            :yhat => OLSout[:yhat],
            :y => Yvec,
            :rsqr => OLSout[:rsqr],
            :rbar => OLSout[:rbar],
            :sige => OLSout[:sige],
            :dw => OLSout[:dw]
        )
    end
    
    # Compute the matrix of coefficients & VCV
    Ft = (X' * X) \ (X' * Y)
    VAR[:Ft] = Ft
    VAR[:F] = Ft'
    SIGMA = (1/(nobse-ntotcoeff)) * (Y - X*Ft)' * (Y - X*Ft)
    VAR[:sigma] = SIGMA
    VAR[:resid] = Y - X*Ft
    VAR[:X] = X
    VAR[:Y] = Y
    if nvar_ex > 0
        VAR[:X_EX] = X_EX
    end
    
    # Companion matrix of F and max eigenvalue
    Fcomp = [VAR[:F][:, 1+constant:nvar*nlag+constant]; 
         [Matrix(I, nvar*(nlag-1), nvar*(nlag-1)) zeros(nvar*(nlag-1), nvar)]]
    VAR[:Fcomp] = Fcomp
    VAR[:maxEig] = maximum(abs.(eigen(Fcomp).values))
    
    # Initialize other results
    VAR[:B] = []
    VAR[:Biv] = []  
    VAR[:PSI] = []
    VAR[:Fp] = []
    VAR[:IV] = []
    
    return VAR, VARopt
end