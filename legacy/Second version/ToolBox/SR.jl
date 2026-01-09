using StatsBase

function SR(VAR, SIGN, VARopt)
    # =======================================================================
    # Compute IRs, VDs, and HDs for a VAR model identified with sign restrictions
    # =======================================================================

    # -------------------- Check inputs --------------------
    if SIGN === nothing
        error("You have not provided sign restrictions (SIGN)")
    end
    if VARopt === nothing
        error("You need to provide VAR options (VARopt from VARmodel)")
    end

    # -------------------- Retrieve parameters --------------------
    nvar    = VAR[:nvar]
    nvar_ex = VAR[:nvar_ex]
    nsteps  = VARopt[:nsteps]    
    ndraws  = VARopt[:ndraws]    
    nobs    = VAR[:nobs]
    nlag    = VAR[:nlag]
    pctg    = VARopt[:pctg]      


    # -------------------- Initialize matrices --------------------
    IRall  = fill(NaN, nsteps, nvar, nvar, ndraws)
    VDall  = fill(NaN, nsteps, nvar, nvar, ndraws)
    HDall = Dict(
        :shock  => zeros(nobs+nlag, nvar, nvar, ndraws),
        :init   => zeros(nobs+nlag, nvar, ndraws),
        :const  => zeros(nobs+nlag, nvar, ndraws),
        :trend  => zeros(nobs+nlag, nvar, ndraws),
        :trend2 => zeros(nobs+nlag, nvar, ndraws),
        :endo   => zeros(nobs+nlag, nvar, ndraws),
        :exo    => zeros(nobs+nlag, nvar, nvar_ex, ndraws)
    )
    Ball = fill(NaN, nvar, nvar, ndraws)

    # -------------------- Sign restriction routine --------------------
    jj = 0   # accepted draws
    tt = 0   # total draws
    ww = 1   # index for printing
    VAR_draw = Dict{String, Any}()

    while jj < ndraws
        if tt > VARopt[:sr_draw]
            println("------------------------------------------------------------")
            println("Total number of draws for finding sign restrictions exceeded")
            println("Change the restrictions or increase VARopt.sr_draw")
            println("------------------------------------------------------------")
            error("See details above")
        end

        # Label for current draw
        label = "draw$(jj)"
        VAR_draw[label] = deepcopy(VAR)
        VARopt[:ident] = "sign"

        # Consider model uncertainty
        if VARopt[:sr_mod] == 1
            sigma_draw, Ft_draw, F_draw, Fcomp_draw = VARdrawpost(VAR)
            VAR_draw[label][:Ft] = Ft_draw
            VAR_draw[label][:F]  = F_draw
            VAR_draw[label][:Fcomp] = Fcomp_draw
            VAR_draw[label][:sigma] = sigma_draw
        end

        # Compute rotated B matrix
        B = SignRestrictions(SIGN, VAR_draw[label], VARopt)

        if !isempty(B)
            # Store accepted draw
            jj += 1
            tt += 1
            Ball[:, :, jj] = B
            VAR_draw[label][:B] = B

            # Compute IR, VD, HD
            aux_irf, VAR_draw[label] = VARir(VAR_draw[label], VARopt)
            IRall[:, :, :, jj] = aux_irf

            aux_fevd, VAR_draw[label] = VARvd(VAR_draw[label], VARopt)
            VDall[:, :, :, jj] = aux_fevd 

            aux_hd, VAR_draw[label] = VARhd(VAR_draw[label], VARopt)
            
            HDall[:shock][:, :, :, jj]  = aux_hd[:shock]
            HDall[:init][:, :, jj]      = aux_hd[:init]
            HDall[:const][:, :, jj]     = aux_hd[:constant]
            HDall[:trend][:, :, jj]     = aux_hd[:trend]
            HDall[:trend2][:, :, jj]    = aux_hd[:trend2]
            HDall[:endo][:, :, jj]      = aux_hd[:endo]
            if nvar_ex > 0
                HDall[:exo][:, :, :, jj] = aux_hd[:endo]
            end

            # Display progress
            if jj == VARopt[:mult] * ww
                println("Rotation: $jj / $ndraws")
                ww += 1
            end
        else
            tt += 1
        end
    end
    println("-- Done!\n")

    # -------------------- Store results --------------------
    SRout = Dict()
    SRout[:IRall] = IRall
    SRout[:VDall] = VDall
    SRout[:Ball]  = Ball
    SRout[:HDall] = HDall

    # Median B
    SRout[:Bmed] = mapslices(median, Ball; dims=3)

    # Find B closest to median
    aux = dropdims(sum((Ball .- SRout[:Bmed]).^2, dims=(1,2)), dims=(1,2))
    sel = findall(aux .== minimum(aux))[1]
    SRout[:B] = Ball[:, :, sel]

    # Compute percentiles
    pctg_inf = (100 - pctg) / 2
    pctg_sup = 100 - (100 - pctg) / 2

    # For impulse responses
    SRout[:IRmed] = mapslices(median, IRall; dims=4)
    SRout[:IRinf] = mapslices(x -> quantile(x, pctg_inf/100), IRall; dims=4)
    SRout[:IRsup] = mapslices(x -> quantile(x, pctg_sup/100), IRall; dims=4)

    # For variance decomposition
    SRout[:VDmed] = mapslices(median, VDall; dims=4)
    SRout[:VDinf] = mapslices(x -> quantile(x, pctg_inf/100), VDall; dims=4)
    SRout[:VDsup] = mapslices(x -> quantile(x, pctg_sup/100), VDall; dims=4)

    # Compute IR, VD, HD for VAR closest to median
    VARopt[:ident] = "sign"
    VAR_sel = VAR_draw["draw$(sel)"]
    SRout[:IR] = VARir(VAR_sel, VARopt)
    SRout[:VD] = VARvd(VAR_sel, VARopt)
    SRout[:HD] = VARhd(VAR_sel, VARopt)

    return SRout
end
