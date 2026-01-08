function VARmakexy(DATA::Matrix{Float64}, lags::Int, constant::Int)
    # =======================================================================
    # Builds VAR matrices from DATA, e.g. [x y] = [x(-1) y(-1) x(-2) y(-2)] 
    # in case of 2 lags and no constant
    # =======================================================================
    
    nobs, nvar = size(DATA)
    
    # Y matrix 
    Y = DATA[lags+1:end, :]
    
    # X-matrix 
    if constant == 0
        X = nothing
        for jj in 0:lags-1
            if X === nothing
                X = DATA[jj+1:nobs-lags+jj, :]
            else
                X = [DATA[jj+1:nobs-lags+jj, :] X]
            end
        end
        
    elseif constant == 1 # constant
        X = nothing
        for jj in 0:lags-1
            if X === nothing
                X = DATA[jj+1:nobs-lags+jj, :]
            else
                X = [DATA[jj+1:nobs-lags+jj, :] X]
            end
        end
        X = [ones(nobs-lags, 1) X]
       
    elseif constant == 2 # time trend and constant
        X = nothing
        for jj in 0:lags-1
            if X === nothing
                X = DATA[jj+1:nobs-lags+jj, :]
            else
                X = [DATA[jj+1:nobs-lags+jj, :] X]
            end
        end
        trend = reshape(collect(1:size(X, 1)), :, 1)  
        X = [ones(nobs-lags, 1) trend X]
     
    elseif constant == 3 # linear time trend, squared time trend, and constant
        X = nothing
        for jj in 0:lags-1
            if X === nothing
                X = DATA[jj+1:nobs-lags+jj, :]
            else
                X = [DATA[jj+1:nobs-lags+jj, :] X]
            end
        end
        trend = reshape(collect(1:size(X, 1)), :, 1)  
        X = [ones(nobs-lags, 1) trend trend.^2 X]
    end
    
    return Y, X
end