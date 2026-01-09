function OLSmodel(y::Vector{Float64}, x::Matrix{Float64}, constant::Int=0)
    # Check inputs
    if isempty(x)
        nobs = length(y)
        nvar = 0
        x = Matrix{Float64}(undef, nobs, 0)  
    else
        nobs, nvar = size(x)
        nobs2 = length(y)
        if nobs != nobs2
            error("x and y must have same # obs")
        end
    end
    
    # Add constant or trend if needed
    if constant == 1 # constant
        x = [ones(nobs, 1) x]
        nvar = nvar + 1
    elseif constant == 2 # trend and constant
        trend = collect(1:nobs)
        x = [ones(nobs, 1) reshape(trend, :, 1) x]
        nvar = nvar + 2
    elseif constant == 3 # trend^2, and constant
        trend = collect(1:nobs)
        x = [ones(nobs, 1) reshape(trend.^2, :, 1) x]
        nvar = nvar + 3
    end
    
    OLS = Dict()
    OLS[:meth] = "ols"
    OLS[:y] = y
    OLS[:x] = x
    OLS[:nobs] = nobs
    OLS[:nvar] = nvar
    
    # xpxi = (X'X)^(-1)
    if nobs < 10000
        Q, R = qr(x)
        xpxi = inv(R' * R)
    else
        xpxi = inv(x' * x)
    end
    
    # OLS estimator
    OLS[:beta] = xpxi * (x' * y)
    
    # Predicted values & residuals
    OLS[:yhat] = x * OLS[:beta]
    OLS[:resid] = y - OLS[:yhat]
    
    # Covariance matrix of residuals
    sigu = OLS[:resid]' * OLS[:resid]
    OLS[:sige] = sigu / (nobs - nvar)
    
    # Covariance matrix of beta
    OLS[:sigbeta] = OLS[:sige] * xpxi
    
    # Std errors of beta, t-stats, intervals, and p-values
    tmp = OLS[:sige] * diag(xpxi)
    sigb = sqrt.(tmp)
    OLS[:bstd] = sigb
    tcrit = -tdis_inv(0.025, Float64(nobs))
    OLS[:bint] = [OLS[:beta] - tcrit .* sigb OLS[:beta] + tcrit .* sigb]
    OLS[:tstat] = OLS[:beta] ./ sqrt.(tmp)
    OLS[:tprob] = tdis_prb(OLS[:tstat], Float64(nobs))
    
    # R2
    ym = y .- mean(y)
    rsqr1 = sigu
    rsqr2 = ym' * ym
    OLS[:rsqr] = 1.0 - rsqr1 / rsqr2
    rsqr1 = rsqr1 / (nobs - nvar)
    rsqr2 = rsqr2 / (nobs - 1.0)
    if rsqr2 != 0
        OLS[:rbar] = 1 - (rsqr1 / rsqr2)
    else
        OLS[:rbar] = OLS[:rsqr]
    end
    
    # Durbin-Watson
    ediff = OLS[:resid][2:end] - OLS[:resid][1:end-1]
    OLS[:dw] = (ediff' * ediff) / sigu
    OLS[:constant] = constant
    
    # F-test
    if constant > 0
        fx = x[:, 1:1] 
        fxpxi = inv(fx' * fx)
        fbeta = fxpxi * (fx' * y)
        fyhat = fx * fbeta
        fresid = y - fyhat
        fsigu = fresid' * fresid
        fym = y .- mean(y)
        frsqr1 = fsigu
        frsqr2 = fym' * fym
        frsqr = 1.0 - frsqr1 / frsqr2
        OLS[:F] = ((frsqr - OLS[:rsqr]) / (nvar - 1)) / ((1 - OLS[:rsqr]) / (nobs - nvar))
    end
    
    return OLS
end