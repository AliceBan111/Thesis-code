function qr_rotation_matrix(K::Int)
    """
    Generate a random K×K orthogonal matrix (equivalent to MATLAB getqr)
    """
    # Draw random matrix
    W = randn(K, K)
    
    # QR decomposition
    F = qr(W)
    Q = Matrix(F.Q)  # Convert Q to regular matrix

    # Ensure diagonal of R is positive
    R = F.R
    for i in 1:K
        if R[i, i] < 0
            Q[:, i] = -Q[:, i]
        end
    end

    return Q
end
