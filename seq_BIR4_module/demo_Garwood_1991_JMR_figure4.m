% demo_Garwood_1991_JMR_figure4.m
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

%% Define parameters for a BIR-4 VS preparation
b1_max = 20;          % maximum RF amplitude [uT]
T_seg  = 2.5e-3;        % duration of one pulse segment [sec]
zeta   = 10;          % constant in the amplitude function [rad]
kappa  = atan(20);    % constant in the frequency/phase function [rad]
beta   = 90;          % flip angle [degree]

%% Calculate dw_max [Hz]
Tp     = 4 * T_seg; % duration of a BIR-4 pulse [sec]
dw_max = 45 * (2 * pi) / Tp;   % maximum frequency sweep [Hz]

%% Calculate a BIR-4 VS module (BIR-4)
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

%% Calculate w1_0 [Hz] and b1_0 [uT]
% b1_0 * gamma * sum(abs(rf_waveform) / b1_max)) * dt = beta
% b1_0 = beta / (gamma * sum(abs(rf_waveform) / b1_max)) * dt)
% [degree] * [cycle/180 degree] / ([Hz/T] * [sec] * [T/1e6uT]) => [uT]
b1_0 = (beta / 180) / (sys.gamma * sum(abs(rf_waveform) / b1_max) * sys.rfRasterTime * 1e-6); % [uT]

% [uT] * [T/1e6uT] * [Hz/T] => *1e-6 [Hz]
w1_0 = (b1_0 * 1e-6) * sys.gamma; % [Hz]

%% Calculate a gradient waveform [mT/m]
grad_waveform = rf_waveform * 0;

%% Perform Bloch simulation
T1 = inf; % T1 relaxation time [sec]
T2 = inf; % T2 relaxation time [sec]

%--------------------------------------------------------------------------
% Calculate the initial magnetization
%--------------------------------------------------------------------------
mx0 = 0;
my0 = 0;
mz0 = 1;

%--------------------------------------------------------------------------
% Calculate the range of off-resonance [Hz]
%--------------------------------------------------------------------------
df_range = (-16:0.1:16).' * w1_0; % [Hz]
Nf = length(df_range);

%--------------------------------------------------------------------------
% Calculate the range of B1 [uT]
%--------------------------------------------------------------------------
B1_range = 1; %(0:0.1:60).' * b1_0; % [uT]
N = length(B1_range);

%--------------------------------------------------------------------------
% Bloch simulation
%--------------------------------------------------------------------------
Mx = zeros(Nf, N, 'double');
My = zeros(Nf, N, 'double');
Mz = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    scale_factor = B1_range(idx) / b1_max;
    [mx,my,mz] = bloch(scale_factor * rf_waveform * 1e-2, grad_waveform * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx(:,idx) = mx;
    My(:,idx) = my;
    Mz(:,idx) = mz;
end

%%

[df_grid, B1_grid] = ndgrid(df_range / w1_0, B1_range / b1_0);
FontSize = 12;

figure('Color', 'w');
surf(df_grid, B1_grid, Mz, 'EdgeColor', 'none');
view(0,-90);
colormap(gray(256));
colorbar;
title('BIR-4 pulse');

%%
Mxy = Mx + 1j * My;
I = double(abs(Mxy) > 0.95);
[df_grid, B1_grid] = ndgrid(df_range / w1_0, B1_range / b1_0);

FontSize = 12;
figure('Color', 'w');
surf(df_grid, B1_grid, I, 'EdgeColor', 'none');
axis tight;
view(0,90);
colormap(parula(256));
colorbar;
title('BIR-4 pulse');
return


%% Display Figure 1
FontSize = 16;

start_index = round(rf_bir4.deadTime / sys.rfRasterTime) + 1;
end_index = find(abs(rf_waveform(start_index:end)) > 0, 1, 'last') + start_index - 1;
index_range = (start_index:end_index).';
N = length(index_range);
t = ((0:N-1).' + 0.5) * sys.rfRasterTime;

figure('Color', 'w', 'Position', [-1 2 828 988]);
color_order = get(gca, 'colororder');
ax1 = subplot(3,1,1);
plot(t * 1e3, abs(rf_waveform(index_range)) * 1e-2, 'Color', color_order(2,:), 'LineWidth', 2);
set(gca, 'FontSize', FontSize, 'XTick' ,(0:2:14).', 'XTickLabel', []);
ylabel('RF (G)', 'FontSize', 16);
xlim([0 t(end)] * 1e3);

ax2 = subplot(3,1,2);
plot(t * 1e3, angle(rf_waveform(index_range)), 'Color', [65 150 81] / 255, 'LineWidth', 2);
set(gca, 'FontSize', FontSize, 'XTick' ,(0:2:14).', 'XTickLabel', [], 'YTickLabel', {'-\pi', '0', '\pi'}, 'YTick', [-pi 0 pi]);
ylabel('Phase (rad)', 'FontSize', 16);
xlim([0 t(end)] * 1e3);
ylim([-5 5]);

% ax3 = subplot(3,1,3);
% plot(t * 1e3, gz_waveform(index_range) * 1e-1, 'Color', color_order(1,:), 'LineWidth', 2);
% set(gca, 'FontSize', FontSize, 'XTick' ,(0:2:14).');
% xlabel('Time (ms)');
% ylabel('G (G/cm)', 'FontSize', 16);
% xlim([0 t(end)] * 1e3);
% ylim([0 4]);
% 
% set(ax1, 'Position', [0.1300 0.7896+0.065*0 0.7750 0.1354]);
% set(ax2, 'Position', [0.1300 0.5705+0.065*1 0.7750 0.1354]);
% set(ax3, 'Position', [0.1300 0.3514+0.065*2 0.7750 0.1354]);
% export_fig(sprintf('Wong_2010_ISMRM_figure1'), '-r300', '-tif');

%% Display Figure 2
[v_grid, B1_grid] = ndgrid(v_range, B1_fraction_range);
FontSize = 12;

figure('Color', 'w', 'Position', [8 590 568 395]);
surf(v_grid, B1_grid * 1e-2, Mz, 'EdgeColor', 'none');
set(gca, 'FontSize', FontSize);
view(0,90);
axis tight;
colormap(jet(256));
hc = colorbar;
title(hc, 'Mz');
xlabel('Velocity (cm/s)');
ylabel('B1 (G)');
export_fig(sprintf('Wong_2010_ISMRM_figure2'), '-r300', '-tif');

%% Display Figure 3
idx1 = find(B1_fraction_range == 0.15 * 1e2);
idx2 = find(B1_fraction_range == 0.20 * 1e2);
idx3 = find(B1_fraction_range == 0.25 * 1e2);

FontSize = 12;
figure('Color', 'w', 'Position', [556 496 684 482]);
color_order = get(gca, 'colororder');
hold on;
plot(v_range, Mz(:,idx1), '-', 'LineWidth', 1, 'Color', color_order(1,:));
plot(v_range, Mz(:,idx2), '--' , 'LineWidth', 1, 'Color', color_order(2,:));
set(gca, 'Box', 'On', 'FontSize', FontSize);
grid on; grid minor;
xlabel('Velocity (cm/s)', 'FontSize', FontSize);
ylabel('Mz', 'FontSize', FontSize);
legend('BIR-4 0.15G', 'BIR-4 0.20G');
export_fig(sprintf('Wong_2010_ISMRM_figure3'), '-r300', '-tif');

return
%% Display a Figure similar to Figure 1
figure('Color', 'w', 'Position', [5 244 1236 734]);
subplot(2,1,1);
yyaxis left;
hold on;
plot(t_RRT * 1e3, real(rf_waveform), 'LineWidth', 1);
plot(t_RRT * 1e3, imag(rf_waveform), 'LineWidth', 1);
set(gca, 'Box', 'On');
ylabel('RF [\muT]', 'FontSize', FontSize);
ylim([-25 25]);

yyaxis right;
plot(t_RRT * 1e3, gz_waveform, 'LineWidth', 1);
ylabel('Gradient amplitude [mT/m]', 'FontSize', FontSize);
xlabel('Time [msec]');
axis tight;
ylim([-50 50]);
grid on; grid minor;
legend('B1+ real', 'B1+ imag', 'Gradient');

subplot(2,1,2);
yyaxis left;
hold on;
plot(t_RRT * 1e3, abs(rf_waveform), 'LineWidth', 1);
set(gca, 'Box', 'On');
ylabel('RF [\muT]', 'FontSize', FontSize);
ylim([-25 25]);

yyaxis right;
plot(t_RRT * 1e3, gz_waveform, 'LineWidth', 1);
ylabel('Gradient amplitude [mT/m]', 'FontSize', FontSize);
xlabel('Time [msec]');
axis tight;
ylim([-50 50]);
grid on; grid minor;
legend('|B1+|', 'Gradient');
