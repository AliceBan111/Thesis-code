function VARmakelags(DATA::Matrix{Float64}, lag::Int)
    # =======================================================================
    # Builds a matrix with lagged values of DATA, i.e. if DATA = [x1 x2],
    # VARmakelags(DATA,1) yields OUT = [x1 x2 x1(-1) x2(-1)]
    # =======================================================================
    
    nobs, nvar = size(DATA)
    
    # Create the lagged matrix
    OUT = Matrix{Float64}(undef, 0, 0)  
    
    for jj in 0:lag-1
        lag_data = DATA[jj+1:nobs-lag+jj, :]
        OUT = [lag_data OUT]  
    end
    
    # Finally, save the non-lagged values...
    aux = DATA[lag+1:end, :]
    
    # ... and append to the lagged matrix
    OUT = [aux OUT]
    
    return OUT
end