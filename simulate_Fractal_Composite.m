% Complete standalone code to generate dataFractal_Composite
% Scenario: Fractal atmospheric phase + localized dense fringes + realistic coherence (geometric decorrelation)
% All required functions are included and can be run directly.

clear; close all; clc;

%% Basic settings
N = 256;        % Image size
L = 1;          % Number of looks (single-look complex data)
rng(0);         % Fix random seed for InSAR noise

%% Generate fractal + localized dense fringe scene
% 1. Generate unwrapped and wrapped phase (fractal background + diagonal ridge)
[psi_frac_unwrapped, psi_frac_wrapped] = make_phase_fractal_local_fringes(N, 1.8, 4*pi);
A_frac = make_amplitude_y_ramp(N, 21, 255);

% 2. Compute realistic coherence using unwrapped phase (capture geometric decorrelation)
rho_frac_real = make_realistic_coherence(N, psi_frac_unwrapped);

% 3. Simulate InSAR interferometric pair (single-look complex data)
dataFractal_Composite = simulate_insar_pair(psi_frac_wrapped, A_frac, rho_frac_real, L);

%% Visualization (optional)
plot_scene(dataFractal_Composite, 'Fractal + Local Fringes + Realistic Coherence');

disp('dataFractal_Composite generated successfully.');

%% ======================== Core Functions ========================

function out = simulate_insar_pair(psi_true, A, rho, L)
    % Simulate single- or multi-look InSAR interferometric pair
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

function [psi_unwrapped, psi_wrapped] = make_phase_fractal_local_fringes(N, alpha, amplitudeScale)
    % Fractal atmospheric phase + localized high-density fringes (ridge)
    
    % 1. Fractal background
    rng(42); 
    noise = randn(N, N);
    F = fftshift(fft2(noise));
    [u, v] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    freq_radius = sqrt(u.^2 + v.^2);
    freq_radius(freq_radius == 0) = 1; 
    F_filtered = F ./ (freq_radius .^ alpha);
    psi_raw = real(ifft2(ifftshift(F_filtered)));
    psi_raw = psi_raw - min(psi_raw(:));
    psi_raw = psi_raw / max(psi_raw(:)); 
    psi_frac = psi_raw * amplitudeScale;      
    
    % 2. Diagonal ridge (producing localized dense fringes)
    [x, y] = meshgrid(1:N, 1:N);
    theta = -pi/4;  % Rotate by 45 degrees
    xr = (x - N/2)*cos(theta) - (y - N/2)*sin(theta);
    yr = (x - N/2)*sin(theta) + (y - N/2)*cos(theta);
    
    % A high, narrow ridge along xr and elongated along yr
    ridge = 25 * pi * exp(-(xr.^2)/(2*(N/10)^2)) .* exp(-(yr.^2)/(2*(N/1.5)^2));
    
    % 3. Superposition
    psi_unwrapped = psi_frac + ridge;
    psi_wrapped = wrapToPi_local(psi_unwrapped);
end

function rho = make_realistic_coherence(N, psi_unwrapped)
    % Simulate realistic coherence map including geometric decorrelation
    
    % 1. Surface texture decorrelation (low-frequency patches)
    rng(100);
    tex_noise = imgaussfilt(randn(N, N), 12);
    tex_noise = tex_noise - min(tex_noise(:));
    tex_noise = tex_noise / max(tex_noise(:));
    rho_texture = 0.3 + 0.6 * tex_noise;   % Coherence in [0.3, 0.9]
    
    % 2. Geometric decorrelation (based on gradient of unwrapped phase)
    [gx, gy] = gradient(psi_unwrapped);
    slope = sqrt(gx.^2 + gy.^2);
    rho_geom = exp(-slope / 2.0);
    
    % 3. Thermal noise fluctuation
    rho_thermal = 0.95 + 0.05 * rand(N, N);
    
    % 4. Combine all effects
    rho = rho_texture .* rho_geom .* rho_thermal;
    rho = min(max(rho, 0.05), 0.99);
end

function A = make_amplitude_y_ramp(N, Amin, Amax)
    % Amplitude map with linear gradient along y-direction
    Ay = linspace(Amin, Amax, N).';
    A  = repmat(Ay, 1, N);
end

function w = wrapToPi_local(x)
    % Phase wrapping to the interval [-pi, pi]
    w = mod(x + pi, 2*pi) - pi;
end

function plot_scene(d, nameStr)
    % Simple visualization tool
    figure('Name', nameStr, 'Color', 'w');
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    
    nexttile; imagesc(d.psi_noisy); axis image off; title([nameStr ' noisy phase']); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(d.psi_true);  axis image off; title([nameStr ' clean phase']); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(d.rho);       axis image off; title('coherence \rho'); colorbar; colormap(gca, 'gray');
    
    nexttile; imagesc(d.A);         axis image off; title('amplitude A'); colorbar; colormap(gca, 'gray');
    nexttile; imagesc(angle(d.I));  axis image off; title('angle(I)'); colorbar; colormap(gca, 'jet');
    nexttile; imagesc(abs(d.I));    axis image off; title('|I|'); colorbar; colormap(gca, 'gray');
end