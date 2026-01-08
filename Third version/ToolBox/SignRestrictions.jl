function SignRestrictions(SIGN, VAR, VARopt)
    # =======================================================================
    # Draws orthonormal rotations of the VAR covariance matrix until 
    # satisfies the signs specified in SIGN.
    # =======================================================================

    dy, ds = size(SIGN)

    # Check whether columns of B (and corresponding sigma) are already provided
    if isnothing(VAR[:Biv]) || isempty(VAR[:Biv])
        if dy == ds
            sigma = VAR[:sigma]
            Biv = nothing
        elseif dy > ds
            sigma = VAR[:sigma]
            aux = cholesky(sigma).L
            Biv = aux[:, dy-ds+1]
        else
            error("Matrix SIGN has the wrong dimension")
        end
    else
        sigma = VAR[:sigma_b]
        Biv = VAR[:Biv]
    end

    # Defines which columns of sigma have to be rotated
    whereToStart = 1 + dy - ds
    if Biv !== nothing && size(Biv, 2) != whereToStart - 1
        error("Biv and SIGN must have coherent sizes")
    end

    nanMat = fill(NaN, dy)
    orderIndices = collect(1:dy)

    # ----------------------------------------------------------------------
    # Search for rotations that satisfy SIGN and Biv
    # ----------------------------------------------------------------------
    counter = 1
    flag = false
    B = nothing
    termaa = zeros(dy, dy)

    while true
        # Create starting matrix to be rotated
        if whereToStart > 1 && Biv !== nothing
            C = cholesky(sigma).L
            q = C \ Biv  # Solve C*q = Biv
            for ii in whereToStart:dy
                r = randn(dy)
                q = hcat(q, ((I - q * q') * r) / norm((I - q * q') * r))
            end
            startingMat = C * q
        else
            startingMat = cholesky(sigma).L
        end

        counter += 1
        if counter > VARopt[:sr_rot]
            println("Max number of rotations reached without satisfying sign restrictions")
            flag = true
            break
        end

        termaa = copy(startingMat)
        TermA = 0

        # Generate random orthogonal matrix for rotation
        rotMat = Matrix{Float64}(I, dy, dy)
        rotMat[whereToStart:end, whereToStart:end] = qr_rotation_matrix(dy - whereToStart + 1)

        # Apply rotation
        terma = termaa * rotMat
        termaa = copy(terma)

        # Update VAR.B for IRF computation
        VAR[:B] = termaa

        # Compute IRF if sr_hor>1
        if VARopt[:sr_hor] > 1
            VARopt_temp = copy(VARopt)
            VARopt_temp[:nsteps] = VARopt[:sr_hor]
            aux_irf, VAR = VARir(VAR, VARopt_temp)
            #VARopt[:nsteps] = VARopt[:sr_hor]
            #aux_irf, VAR = VARir(VAR, VARopt)
        end


        # Check sign restrictions
        for ii in 1:ds
            for jj in whereToStart:dy
                if isfinite(terma[1, jj])
                    CHECK = VARopt[:sr_hor] > 1 ? permutedims(aux_irf[:, :, jj]) : terma[:, jj]
                    if sum(CHECK .* SIGN[:, ii] .< 0) == 0
                        TermA += 1
                        orderIndices[whereToStart-1+ii] = jj
                        terma[:, jj] .= nanMat
                        break
                    elseif sum(-CHECK .* SIGN[:, ii] .< 0) == 0
                        TermA += 1
                        terma[:, jj] .= nanMat
                        termaa[:, jj] .= -termaa[:, jj]
                        orderIndices[ii+whereToStart-1] = jj
                        break
                    end
                end
            end
        end
        if TermA == ds
            break
        end
    end

    if !flag
        B = termaa[:, orderIndices]
    else
        B = Array{Float64}(undef, 0, 0)
    end

    return B
end
