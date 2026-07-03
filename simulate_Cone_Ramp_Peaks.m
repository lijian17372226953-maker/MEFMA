% Complete MATLAB code to generate simulated InSAR data 
% Covers 3 Scenarios: Cone / Peaks / Octant Ramp
% clear; close all; clc;

%% Settings
N = 256;        % image size
L = 1;             % number of looks (set >1 for multilook)
rng(0);           % reproducibility

% Coherence map: 0.1 -> 0.9 left-to-right
rho = repmat(linspace(0.1, 0.9, N), N, 1);

%% ===== Scene 1: CONE =====
psi_cone = make_phase_cone(N, 24);
A_cone   = make_amplitude_y_ramp(N, 21, 255);
dataCone = simulate_insar_pair(psi_cone, A_cone, rho, L);

%% ===== Scene 2: PEAKS (MATLAB peaks) =====
psi_peaks = make_phase_peaks_matlab(N, 5);
A_peaks   = make_amplitude_y_ramp(N, 21, 255);
dataPeaks = simulate_insar_pair(psi_peaks, A_peaks, rho, L);

%% ===== Scene 3: OCTANT RAMP =====
psi_ramp = make_phase_octant_ramp(N, 28, 8);
A_ramp   = 128 * ones(N);
dataRamp = simulate_insar_pair(psi_ramp, A_ramp, rho, L);

%% Visualize
plot_scene(dataCone,    '1. Cone');
plot_scene(dataPeaks,   '2. Peaks');
plot_scene(dataRamp,    '3. Octant Ramp (8 Slices)');

disp('Done. Simulated data generated and visualized.');

%% ======================= Core simulator ============================
function out = simulate_insar_pair(psi_true, A, rho, L)
    N = size(psi_true,1);
    psi_true = wrapToPi_local(psi_true);
    rho = min(max(rho, 0), 0.9999); 
    
    psi3 = repmat(psi_true, 1, 1, L);
    A3   = repmat(A,        1, 1, L);
    rho3 = repmat(rho,      1, 1, L);
    
    v1 = (randn(N,N,L) + 1i*randn(N,N,L)) / sqrt(2);
    v2 = (randn(N,N,L) + 1i*randn(N,N,L)) / sqrt(2);
    
    u1 = A3 .* v1;
    u2 = A3 .* ( rho3 .* exp(-1i*psi3) .* v1 + sqrt(1 - rho3.^2) .* v2 );
    
    Ik = u1 .* conj(u2);
    I  = mean(Ik, 3);
    
    out.psi_true  = psi_true;
    out.psi_noisy = angle(I);
    out.rho       = rho;
    out.A         = A;
    out.L         = L;
    out.u1 = u1(:,:,1);
    out.u2 = u2(:,:,1);
    out.I  = I;
end

%% ======================= Phase fields ==============================
function psi = make_phase_cone(N, radialPeriodPx)
    [x,y] = meshgrid(1:N, 1:N);
    cx = (N+1)/2; cy = (N+1)/2;
    r  = sqrt((x-cx).^2 + (y-cy).^2);
    psi = 2*pi * (r / radialPeriodPx);
    psi = wrapToPi_local(psi);
end

function psi = make_phase_peaks_matlab(N, wraps)
    Z = peaks(N);
    Z = Z - mean(Z(:));
    Z = Z / (max(abs(Z(:))) + 1e-12);   
    psi = (wraps * 2*pi) * Z;
    psi = wrapToPi_local(psi);
end

function psi = make_phase_octant_ramp(N, periodBottom, periodTop)
    cx = (N+1)/2;
    cy = (N+1)/2;
    margin = N; 
    y_ext = (1 - margin : N + margin).';
    P_ext = periodTop + (periodBottom - periodTop) * (y_ext - 1) / (N - 1);
    P_ext = max(P_ext, 2); 
    dphi = 2*pi ./ P_ext;
    phi_ext = cumsum(dphi);
    phi_ext = phi_ext - interp1(y_ext, phi_ext, cy);
    [X, Y] = meshgrid(1:N, 1:N);
    Xc = X - cx;
    Yc = Y - cy;
    theta = mod(atan2(Yc, Xc), 2*pi);
    k = floor(theta / (pi/4));
    gamma = k * (pi/4);
    Yrot = -Xc .* sin(gamma) + Yc .* cos(gamma);
    Y_sample = Yrot + cy;
    psi_1d = interp1(y_ext, phi_ext, Y_sample(:), 'linear', 'extrap');
    psi = reshape(psi_1d, N, N);
    psi = wrapToPi_local(psi);
end

%% ======================= Amplitude fields ==========================
function A = make_amplitude_y_ramp(N, Amin, Amax)
    Ay = linspace(Amin, Amax, N).';
    A  = repmat(Ay, 1, N);
end

%% ======================= Utilities ================================
function w = wrapToPi_local(x)
    w = mod(x + pi, 2*pi) - pi;
end

function plot_scene(d, nameStr)
    figure('Name', nameStr, 'Color', 'w', 'Position', [100, 100, 1200, 600]);
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    
    nexttile; imagesc(d.psi_noisy); axis image off; title([nameStr ' noisy phase']); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(d.psi_true);  axis image off; title([nameStr ' clean phase']); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(d.rho);       axis image off; title('coherence \rho'); colorbar; colormap(gca, 'gray');
    
    nexttile; imagesc(d.A);         axis image off; title('amplitude A'); colorbar; colormap(gca, 'gray');
    nexttile; imagesc(angle(d.I));  axis image off; title('angle(I)'); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(abs(d.I));    axis image off; title('|I|'); colorbar; colormap(gca, 'gray');
end