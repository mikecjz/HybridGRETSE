% demo_Cunningham_2006_MRM_figure1.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/08/2022, Last modified: 09/08/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
src_directory = '';
thirdparty_directory = 'D:\VSASL\thirdparty';
pulseq_directory = 'D:\pulseq\pulseq';

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(pulseq_directory));

%% Define imaging parameters
% 180.18 / 100/ 50 =  fast / normal / whisper
Gmax = 30;   % max gradient strength [mT/m]
Smax = 40;   % maximum slew rate [mT/m/ms]
B0   = 0.55; % main field strength [T]

%% Set system limits
sys = mr.opts('MaxGrad', Gmax, 'GradUnit', 'mT/m' , ...
              'MaxSlew', Smax, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6 , ...
              'rfDeadtime'    , 100e-6, ...
              'adcDeadTime'   , 10e-6 , ...
              'B0', B0);

%% Define parameters for a BIR-4 preparation
T_seg     = 2e-3;        % duration of one pulse segment [sec]
b1_max    = 20;          % maximum RF amplitude [uT]
dw_max    = 39.8e3;      % maximum frequency sweep [Hz]
zeta      = 15.2;        % constant in the amplitude function [rad]
kappa     = atan(63.6);  % constant in the frequency/phase function [rad]
beta      = 90;          % flip angle [degree]

%--------------------------------------------------------------------------
% Comparison between Garwood 1991 JMR and Staewen 1989 IR
%--------------------------------------------------------------------------
% Garwood 1991 JMR:
% For Segment 1 (0 < t <= 0.25 * Tp)
% w1(t) = w1_max * tanh(zeta * (1 - 4 * t / Tp))
%
% Staewen 1989 IR:
% segment 1: f_B(t) = tanh(beta * (1 - tau)) (0 < tau < 1)
%
% Starting from (0 < t <= 0.25 * Tp)
% (0 < 4 * t <= Tp)
% (0 < 4 * t / Tp <= 1) := (0 < tau <= 1)
%    w1(t) = w1_max * tanh(zeta * (1 - 4 * t / Tp))
% => w1(t) = w1_max * tanh(zeta * (1 - tau))
%
% variables in Garwood 1991 JMR vs variables in Staewen 1989 IR
% zeta == beta
%--------------------------------------------------------------------------
% kappa_staewen = atan(10);
% beta_staewen = 1.6; % From Figure 1 of Cunningham 2006 MRM
% 
% kappa_staewen = atan(63.6);
% beta_staewen = 15.2;
% 
% dw_max_staewen = 100 * pi / (4 * T_seg) / (2 * pi); % [Hz]
% lambda = 5 * pi;
% zeta = beta_staewen; 
% kappa = kappa_staewen;
% dw_max = dw_max_staewen;

%% Calculate a BIR-4 module (BIR-4)
rf_bir4 = calculate_pulseq_BIR4_module(T_seg, b1_max, dw_max, zeta, kappa, beta, sys);

%% Create a sequence object
seq = mr.Sequence(sys);

%% Add a new block to the sequence
seq.addBlock(rf_bir4);

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 2] * 1);

%% Get discrete samples
wave_data = seq.waveforms_and_times(true);

%% Interpolate RF and gradient waveforms to RRT
% [Hz] / [Hz/T] * [1e6uT/T] => *1e6 [uT]
N_RRT = round(mr.calcDuration(rf_bir4) / sys.rfRasterTime);
t_RRT = ((0:N_RRT-1).' + 0.5) * sys.rfRasterTime;

%--------------------------------------------------------------------------
% Interpolate an RF waveform to RRT [uT]
%--------------------------------------------------------------------------
t_rf = cat(1, t_RRT(1), wave_data{4}(1,:).', t_RRT(end));
rf = cat(1, 0, conj(wave_data{4}(2,:)).' / sys.gamma * 1e6, 0); % [uT]
rf_waveform = interp1(t_rf, rf, t_RRT, 'linear', 'extrap');

%% Scale the RF waveform according to lambda = gamma * sum(|B1(t)|) * dt [rad]
% [Hz/T] * [2pi rad/cycle] * [T/1e6uT] * [uT] * [sec] => [rad]
%current_lambda = (sys.gamma * 2 * pi * 1e-6) * sum(abs(rf_waveform)) * sys.rfRasterTime;
%rf_waveform = rf_waveform / current_lambda * lambda;

%% Calculate a gradient waveform [mT/m]
grad_waveform = rf_waveform * 0;

%% Perform Bloch simulation
T1 = inf; % T1 relaxation time [sec]
T2 = inf; % T2 relaxation time [sec]

%--------------------------------------------------------------------------
% Calculate the range of off-resonance [Hz]
%--------------------------------------------------------------------------
Nf = 400;
df_range = linspace(-500, 500, Nf); % [Hz]

%--------------------------------------------------------------------------
% Calculate the range of B1 scales
%--------------------------------------------------------------------------
N = 400;
B1_scale_range = linspace(0, 2, N);

%--------------------------------------------------------------------------
% Bloch simulation for Figure 1d
%--------------------------------------------------------------------------
mx0 = zeros(1, Nf, 'double');
my0 = zeros(1, Nf, 'double');
mz0 = ones(1, Nf, 'double');

Mx1 = zeros(Nf, N, 'double');
My1 = zeros(Nf, N, 'double');
Mz1 = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    B1_scale = B1_scale_range(idx);
    [mx,my,mz] = bloch(B1_scale * rf_waveform * 1e-2, grad_waveform * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx1(:,idx) = mx;
    My1(:,idx) = my;
    Mz1(:,idx) = mz;
end

%% Bloch simulation for Figure 1e
mx0 = ones(1, Nf, 'double');
my0 = zeros(1, Nf, 'double');
mz0 = zeros(1, Nf, 'double');

Mx2 = zeros(Nf, N, 'double');
My2 = zeros(Nf, N, 'double');
Mz2 = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    B1_scale = B1_scale_range(idx);
    [mx,my,mz] = bloch(B1_scale * rf_waveform * 1e-2, grad_waveform * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx2(:,idx) = mx;
    My2(:,idx) = my;
    Mz2(:,idx) = mz;
end

%% Bloch simulation for Figure 1f
mx0 = zeros(1, Nf, 'double');
my0 = ones(1, Nf, 'double');
mz0 = zeros(1, Nf, 'double');

Mx3 = zeros(Nf, N, 'double');
My3 = zeros(Nf, N, 'double');
Mz3 = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    B1_scale = B1_scale_range(idx);
    [mx,my,mz] = bloch(B1_scale * rf_waveform * 1e-2, grad_waveform * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx3(:,idx) = mx;
    My3(:,idx) = my;
    Mz3(:,idx) = mz;
end

%% Display Figure 1(d,e,f)
[df_grid, B1_grid] = ndgrid(df_range, B1_scale_range);

FontSize = 12;
figure('Color', 'w', 'Position', [402 -1 425 971]);
ax1 = subplot(3,1,1);
surf(df_grid, B1_grid, abs(Mz1), 'EdgeColor', 'none');
view(0,-90);
axis square;
colormap(gray(256));
hc1 = colorbar;
caxis([0 1]);
title('BIR-4 pulse');
text(-440, 1.8, 'M(0) = (0,0,M_0)', 'Color', 'w', 'FontSize', FontSize, 'FontAngle', 'italic');

ax2 = subplot(3,1,2);
surf(df_grid, B1_grid, abs(Mz2), 'EdgeColor', 'none');
view(0,-90);
axis square;
colormap(gray(256));
hc2 = colorbar;
caxis([0 1]);
text(-440, 1.8, 'M(0) = (M_0,0,0)', 'Color', 'w', 'FontSize', FontSize, 'FontAngle', 'italic');

ax3 = subplot(3,1,3);
surf(df_grid, B1_grid, abs(Mz3), 'EdgeColor', 'none');
view(0,-90);
axis square;
colormap(gray(256));
hc3 = colorbar;
caxis([0 1]);
text(-440, 1.8, 'M(0) = (0,M_0,0)', 'Color', 'k', 'FontSize', FontSize, 'FontAngle', 'italic');

set(ax1, 'Position', [0.1300 0.7093 - 0.04 0.5742 + 0.06 0.2157 + 0.06]);
set(ax2, 'Position', [0.1300 0.4096 - 0.06 0.5742 + 0.06 0.2157 + 0.06]);
set(ax3, 'Position', [0.1300 0.1100 - 0.08 0.5742 + 0.06 0.2157 + 0.06]);

set(hc1, 'Position', [0.8501 - 0.07 0.6693 0.0502 0.2757]);
set(hc2, 'Position', [0.8501 - 0.07 0.3496 0.0502 0.2757]);
set(hc3, 'Position', [0.8501 - 0.07 0.0300 0.0502 0.2757]);

export_fig(sprintf('Cunningham_2006_MRM_figure1'), '-r300', '-tif');
