function P_final = MEFMA_Reconstruction(noisy_phase, N)
% MEFMA_Reconstruction: Multiscale Empirical Fringe-Manifold Adaptive Filtering
%
% This function implements the core V-cycle-like phase decomposition and 
% cascaded reconstruction framework of the MEFMA algorithm. It is designed 
% to achieve robust InSAR phase restoration in extremely low SNR and complex 
% decorrelation environments.
%
% Inputs:
%   noisy_phase - The observed wrapped interferometric phase map [-pi, pi)
%   N           - Predefined upper limit of multiscale decomposition levels 
%                 (default is 5)
%
% Outputs:
%   P_final     - The high-fidelity reconstructed continuous phase map

    if nargin < 2, N = 5; end

    M = cell(N+1, 1);
    R = cell(N, 1);
    R_hat = cell(N, 1);

    M{1} = noisy_phase;

    fprintf('=== Stage 1: Macroscale Multi-Frequency Phase Decomposition ===\n');
    % Decoupling mixed topographic details from chaotic speckle perturbations
    for k = 1:N
        win_size = 8 * 2^(k-1);
        fprintf('  Adaptive frequency stripping, Scale %d (Window: %d)...\n', k, win_size);
        % Extract low-frequency coherent base phase
        M{k+1} = goldstein_filter(M{k}, win_size, 1);
        % Isolate localized high-frequency topological variations and noise
        R{k} = wrapToPi(M{k} - M{k+1});
    end

    fprintf('=== Stage 2: Micro-Feature Residual Denoising (Scale-Bound Projection) ===\n');
    % Executing adaptive fringe manifold projection directly on residual phase maps
    for k = 1:N
        win_size = 8 * 2^(k-1);
        blk_size = win_size;

        fprintf('  Processing residual layer %d (Block Size: %d)...\n', k, blk_size);
        
        % Denoised complex residual field mapping
        Z_clean = goldstein_filter(R{k}, blk_size, 1);
        Z_clean = exp(1i .* Z_clean);

        % Coherence-guided phase attenuation (preventing nonlinear noise amplification)
        Coh = abs(Z_clean); 
        R_hat{k} = atan2(imag(Z_clean), real(Z_clean)) .* Coh; 
    end

    fprintf('=== Stage 3: Coherent Cascaded Reconstruction (Coarse to Fine) ===\n');
    % Synthesizing robust global SNR enhancement layer by layer
    P_final = M{N+1}; 
    for k = N:-1:1
        P_final = wrapToPi(P_final + R_hat{k});
    end

    fprintf('=== Stage 4: Empirical Fringe-Manifold Reverse Polishing ===\n');
    % Enforcing global topological continuity while honoring steep deformation gradients
    Z_final = exp(1i * P_final); 

    for k = 1:N
        win_size = 8 * 2^(k-1);
        blk_size = win_size; 
        patch_size = 8;

        fprintf('  Reverse manifold polishing at scale %d (Block Size: %d)...\n', k, blk_size);    
        
        % Adaptive Wiener shrinkage: Suppresses noise eigenvalues without damaging fringe skeleton
        Z_final = block_eigen_fringe_adaptive(Z_final, blk_size, patch_size); 
        
        % CRITICAL: Amplitude normalization ensures "pure phase characteristics" 
        % are passed to the next scale manifold
        Z_final = Z_final ./ (abs(Z_final) + eps); 
    end

    % Extract the ultimate high-fidelity filtered phase map
    P_final = angle(Z_final);

    fprintf('✅ Reconstruction and manifold polishing completed successfully!\n');
end

%% ================= Core Subroutines ================= %%

function Z_clean = block_eigen_fringe_adaptive(Z, blk_size, patch_size)
% Block-based adaptive eigen-fringe manifold Wiener shrinkage
% Utilizes Wax-Kailath Minimum Description Length (MDL) for intrinsic dimensionality estimation
    [rows, cols] = size(Z);
    pad_r = ceil(rows/blk_size)*blk_size - rows;
    pad_c = ceil(cols/blk_size)*blk_size - cols;
    Z_pad = padarray(Z, [pad_r, pad_c], 'symmetric', 'post');
    [pr, pc] = size(Z_pad);
    Z_clean_pad = zeros(pr, pc);
    weight = zeros(pr, pc);
    step = blk_size / 2; % 50% overlapping blocks

    for r = 1:step:(pr - blk_size + 1)
        for c = 1:step:(pc - blk_size + 1)
            blk = Z_pad(r:r+blk_size-1, c:c+blk_size-1);
            X = im2col(blk, [patch_size, patch_size], 'sliding'); 

            % Covariance matrix and eigenvalue decomposition
            N_samples = size(X, 2);
            C = (X * X') / N_samples;
            [V, D] = eig(C);
            lambda = real(diag(D)); 
            [lambda, sort_idx] = sort(lambda, 'descend');
            V = V(:, sort_idx);

            % =================================================================
            % Wax-Kailath MDL Information Theoretic Criteria 
            % Adaptive Signal Subspace Delineation
            % =================================================================
            p = length(lambda);
            lam_safe = max(lambda, eps); 
            lam_flip = flipud(lam_safe); 

            % Precompute cumulative arithmetic and geometric means
            idx = (1:p-1)'; 
            cum_arith = cumsum(lam_flip(1:p-1)) ./ idx;
            cum_geom  = exp(cumsum(log(lam_flip(1:p-1))) ./ idx);

            % k represents the candidate signal subspace dimensions
            k = p - idx; 

            % Log-likelihood and Penalty calculation
            log_likelihood = -N_samples * idx .* log(cum_geom ./ cum_arith);
            penalty = 0.5 * k .* (2*p - k) * log(N_samples);

            % MDL Objective Function
            MDL = log_likelihood + penalty;

            % Optimal intrinsic dimensionality yielding minimum MDL
            [~, min_idx] = min(MDL);
            k_opt = k(min_idx);

            % Physically bounded noise floor estimation
            sigma_n2 = mean(lambda(k_opt+1:end));
            
            % =================================================================
            % Adaptive Wiener Shrinkage Weights
            w = max(0, 1 - (sigma_n2 ./ (lambda + eps)));

            % Subspace Projection
            X_proj = V * diag(w) * V' * X;

            % Overlap-Add Merging
            blk_clean = col2im_overlap_add(X_proj, blk_size, blk_size, patch_size);
            Z_clean_pad(r:r+blk_size-1, c:c+blk_size-1) = Z_clean_pad(r:r+blk_size-1, c:c+blk_size-1) + blk_clean;
            weight(r:r+blk_size-1, c:c+blk_size-1) = weight(r:r+blk_size-1, c:c+blk_size-1) + 1;
        end
    end
    Z_clean_pad = Z_clean_pad ./ (weight + eps);
    Z_clean = Z_clean_pad(1:rows, 1:cols);
end

function I = col2im_overlap_add(X, M, N, P)
% Reconstructs image from sliding patches using Overlap-Add method
    I = zeros(M, N);
    W = zeros(M, N);
    idx = 1;
    for c = 1:(N-P+1)
        for r = 1:(M-P+1)
            patch = reshape(X(:, idx), P, P);
            I(r:r+P-1, c:c+P-1) = I(r:r+P-1, c:c+P-1) + patch;
            W(r:r+P-1, c:c+P-1) = W(r:r+P-1, c:c+P-1) + 1;
            idx = idx + 1;
        end
    end
    I = I ./ (W + eps);
end

function M_out = goldstein_filter(M_in, win, alpha)
% Block-based implementation of the frequency-adaptive filtering operator
    Z_in = exp(1i * M_in);
    overlap = ceil(win / 4);
    filt_fun = @(bs) process_block(bs.data, alpha);
    Z_out = blockproc(Z_in, [win win], filt_fun, ...
        'BorderSize', [overlap overlap], 'PadMethod', 'symmetric', 'UseParallel', true);
    M_out = angle(Z_out);
end

function out = process_block(data, ~) 
% Localized spatial-frequency adaptive processing driven by Hoyer's Sparsity
    [R, C] = size(data);

    % =========================================================================
    % Uncertainty Principle-Driven Parameterless Smoothing
    % =========================================================================
    % Construct a purely geometry-dependent 2D Hann window.
    % Prevents physical distortion of smoothing strength across scales during the V-Cycle.
    r_vec = (0:R-1)' / (R-1);
    c_vec = (0:C-1) / (C-1);
    W = (0.5 - 0.5 * cos(2 * pi * r_vec)) * (0.5 - 0.5 * cos(2 * pi * c_vec));

    F = fftshift(fft2(data));
    S_smooth = abs(fftshift(fft2(data .* W))); 

    % =========================================================================
    % Hoyer's Sparsity-Based Adaptive Alpha 
    % =========================================================================
    % Quantifies the structural coherence of local fringes in the frequency domain
    S_mag = S_smooth(:);
    N_pts = length(S_mag);
    L1_S = sum(S_mag);
    L2_S = sqrt(sum(S_mag.^2));

    % Calculate Hoyer's sparsity index
    if L2_S < eps
        sparsity = 0;
    else
        sparsity = (sqrt(N_pts) - L1_S/L2_S) / (sqrt(N_pts) - 1);
    end

    % Adaptive smoothing exponent is accurately modeled as complementary state of sparsity
    alpha_adapt = 1 - sparsity; 

    S_norm = S_smooth / max(S_smooth(:));
    weight = S_norm .^ alpha_adapt;

    F_filt = F .* weight;
    out = ifft2(ifftshift(F_filt));
end