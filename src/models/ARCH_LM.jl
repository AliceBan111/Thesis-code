using HypothesisTests, GLM, DataFrames, Distributions

"""
对VAR残差的每一列做ARCH-LM检验
resids: T×K 残差矩阵
lags:   检验的滞后阶数
"""
function ARCH_LM_test(resids::Matrix{Float64}; lags::Int=4)
    T, K = size(resids)
    println("=== ARCH-LM Test (lags=$lags) ===")
    
    for k in 1:K
        e = resids[:, k]
        e2 = e .^ 2
        
        # 构造回归：e²_t = a0 + a1*e²_{t-1} + ... + ap*e²_{t-p}
        Y = e2[lags+1:end]
        X = hcat([e2[lags+1-j:end-j] for j in 1:lags]...)
        X = hcat(ones(length(Y)), X)
        
        # OLS
        b = X \ Y
        Yhat = X * b
        
        # LM统计量 = T * R²
        SS_res = sum((Y .- Yhat).^2)
        SS_tot = sum((Y .- mean(Y)).^2)
        R2 = 1 - SS_res / SS_tot
        LM = length(Y) * R2
        
        # χ²(lags)分布
        pval = 1 - cdf(Chisq(lags), LM)
        
        println("Variable $k: LM = $(round(LM, digits=4)), p-value = $(round(pval, digits=4))")
        pval < 0.05 ? println("reject H0, ARCH effect") : println("not reject H0, no ARCH effect")
    end
end

function plot_resid_diagnostics(resids::Matrix{Float64}, varnames::Vector{String})
    K = size(resids, 2)
    plt = plot(layout=(K, 2), size=(1000, 200*K))

    for k in 1:K
        e = resids[:, k]
        # 左：残差时序
        plot!(plt, subplot=2k-1, e, 
              title="$(varnames[k]) resid", 
              legend=false)
        hline!(plt, [0.0], subplot=2k-1, linestyle=:dash, color=:black)
        
        # 右：残差平方（看波动聚集）
        plot!(plt, subplot=2k, e.^2, 
              title="$(varnames[k]) resid²", 
              legend=false, color=:orange)
    end
    display(plt)
end

