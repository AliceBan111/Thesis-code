function CommonSample(DATA::Matrix{Float64}, dim::Int=1)
    # =======================================================================
    # If a row of DATA contains a NaN, the row is removed. If dim=2, if a 
    # column of DATA contains a NaN, the column is removed.
    # =======================================================================
    
    fo = 0
    lo = 0

    if dim == 1
        temp = vec(sum(DATA, dims=2))  
        ii = 1
        if isnan(temp[ii])
            while isnan(temp[ii])
                fo = fo + 1
                ii = ii + 1
                if ii > length(temp)
                    break
                end
            end
        end
        for ii in 1:size(DATA,1)-fo
            if isnan(temp[end+1-ii])
                lo = lo + 1
            end
        end
        DATA = DATA[vec(.!any(isnan.(DATA), dims=2)), :]
    else
        temp = vec(sum(DATA, dims=1)) 
        ii = 1
        if isnan(temp[ii])
            while isnan(temp[ii])
                fo = fo + 1
                ii = ii + 1
            end
        end
        for ii in 1:size(DATA,2)-fo
            if isnan(temp[end+1-ii])
                lo = lo + 1
            end
        end
        DATA = DATA[:, vec(.!any(isnan.(DATA), dims=1))]
    end

    return DATA, fo, lo
end