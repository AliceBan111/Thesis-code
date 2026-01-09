"""
    extract_structural_shocks(res_dict; shock_index=1)

Extracts structural shocks from a VAR result dictionary.
- `res_dict`: The result dictionary returned by the `estimate_VAR_SR` function.
- `shock_index`: The column index of the target structural shock (defaults to 1, typically the Price Markup Shock).
"""
function extract_structural_shocks(res_dict; shock_index=1)
    # 1. Integrity check: ensure required sub-dictionaries exist
    if !haskey(res_dict, :VAR) || !haskey(res_dict, :SRout)
        error("KeyError: The result dictionary must contain :VAR and :SRout keys.")
    end
    
    var_part = res_dict[:VAR]
    sr_part = res_dict[:SRout]
    
    # 2. Extract the reduced-form residual matrix (Size: (T-p) x N)
    resid_matrix = var_part[:resid]
    
    # 3. Retrieve the median impact matrix B from accepted draws (Size: N x N)
    ball = sr_part[:Ball]
    n_accepted = size(ball, 3)
    median_idx = div(n_accepted, 2)
    B_median = ball[:, :, median_idx]
    
    # 4. Compute structural shocks
    # Formula: Structural Shocks = Residuals * inv(B_median)'
    # This transforms reduced-form errors into orthogonal structural shocks.
    all_structural_shocks = resid_matrix * inv(B_median)'
    
    # 5. Return the specific shock series
    return all_structural_shocks[:, shock_index]
end