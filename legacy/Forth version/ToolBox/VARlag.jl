function VARlag(ENDO::Matrix{Float64}, maxlag::Int, constant::Int=1, EXOG::Union{Matrix{Float64}, Nothing}=nothing, lag_ex::Int=0)
    # =======================================================================
    # Determine VAR lag length with Akaike (AIC) and Schwarz Bayesian 
    # Criterion (SBC) criterion.
    # =======================================================================
    
    ## Check inputs
    nobs, _ = size(ENDO)
    
    # Check if there are exogenous variables
    if EXOG !== nothing
        nobs2, num_ex = size(EXOG)
        # Check that ENDO and EXOG are conformable
        if nobs2 != nobs
            error("var: nobs in EXOG-matrix not the same as y-matrix")
        end
    else
        num_ex = 0
    end
    
    # number of exogenous variables per equation
    nvar_ex = num_ex * (lag_ex + 1)
    
    ## Compute log likelihood and Akaike criterion
    logL = zeros(maxlag)
    AIC_vec = zeros(maxlag)
    SBC_vec = zeros(maxlag)
    
    for i in 1:maxlag
        X = ENDO[maxlag+1-i:end, :]

        aux, _ = VARmodel(X, i, constant)

        if nvar_ex > 0
            Y = EXOG[maxlag+1-i:end, :]
            aux, _ = VARmodel(X, i, constant, Y, lag_ex)
        end
        
        NOBSadj = aux[:nobs]
        NOBS = aux[:nobs] + i  
        NVAR = aux[:nvar]
        NTOTCOEFF = aux[:ntotcoeff]
        RES = aux[:resid]
        
        # VCV of the residuals (use dof adjusted denominator)
        SIGMA = (1/NOBSadj) * (RES' * RES)
        
        # Log-likelihood
        logL[i] = -(NOBS/2) * (NVAR*(1 + log(2*π)) + log(det(SIGMA)))
        
        # AIC: −2*LogL/T + 2*n/T, where n is total number of parameters (ie, NVAR*NTOTCOEFF)
        AIC_vec[i] = -2*(logL[i]/NOBS) + 2*(NVAR*NTOTCOEFF)/NOBS
        
        # SBC: −2*LogL/T + n*log(T)/T
        SBC_vec[i] = -2*(logL[i]/NOBS) + (NVAR*NTOTCOEFF)*log(NOBS)/NOBS
    end
    
    # Find the min of the info criteria
    AIC = argmin(AIC_vec)
    SBC = argmin(SBC_vec)
    
    return AIC, SBC, logL
end