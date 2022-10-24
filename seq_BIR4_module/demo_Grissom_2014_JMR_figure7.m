% demo_Grissom_2014_JMR_figure7.m
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
zeta   = 10;                       % constant in the amplitude function [rad]
beta   = 45;                       % flip angle [degree]
b1_max = 1;                        % maximum RF amplitude [uT]

% 4.7-ms BIR-4 pulse
Tp1     = 4.7e-3;                    % duration of a BIR-4 pulse [sec]
T_seg1  = Tp1 / 4;                   % duration of one pulse segment [sec]
dw_max1 = 100 * pi / Tp1 / (2 * pi); % maximum frequency sweep [Hz]
kappa1  = atan(20);                  % constant in the frequency/phase function [rad]

% 5.9-ms BIR-4 pulse
Tp2     = 5.9e-3;                    % duration of a BIR-4 pulse [sec]
T_seg2  = Tp2 / 4;                   % duration of one pulse segment [sec]
dw_max2 = 100 * pi / Tp2 / (2 * pi); % maximum frequency sweep [Hz]
kappa2  = atan(15);                  % constant in the frequency/phase function [rad]

%% Calculate a BIR-4 VS module (BIR-4)
rf1_bir4 = calculate_pulseq_BIR4_module(T_seg1, b1_max, dw_max1, zeta, kappa1, beta, sys);
rf2_bir4 = calculate_pulseq_BIR4_module(T_seg2, b1_max, dw_max2, zeta, kappa2, beta, sys);


return

%% Get discrete samples
seq = mr.Sequence(sys);
seq.addBlock(rf1_bir4);
wave_data1 = seq.waveforms_and_times(true);

seq = mr.Sequence(sys);
seq.addBlock(rf2_bir4);
wave_data2 = seq.waveforms_and_times(true);

%% Interpolate RF and gradient waveforms to RRT
% [Hz] / [Hz/T] * [1e6uT/T] => *1e6 [uT]
N_RRT1 = round(mr.calcDuration(rf1_bir4) / sys.rfRasterTime);
N_RRT2 = round(mr.calcDuration(rf2_bir4) / sys.rfRasterTime);
t_RRT1 = ((0:N_RRT1-1).' + 0.5) * sys.rfRasterTime;
t_RRT2 = ((0:N_RRT2-1).' + 0.5) * sys.rfRasterTime;

%--------------------------------------------------------------------------
% Interpolate an RF waveform to RRT [uT]
%--------------------------------------------------------------------------
t_rf1 = cat(1, t_RRT1(1), wave_data1{4}(1,:).', t_RRT1(end));
rf1 = cat(1, 0, conj(wave_data1{4}(2,:)).' / sys.gamma * 1e6, 0); % [uT]
rf_waveform1 = interp1(t_rf1, rf1, t_RRT1, 'linear', 'extrap');

%--------------------------------------------------------------------------
% Interpolate an RF waveform to RRT [uT]
%--------------------------------------------------------------------------
t_rf2 = cat(1, t_RRT2(1), wave_data2{4}(1,:).', t_RRT2(end));
rf2 = cat(1, 0, conj(wave_data2{4}(2,:)).' / sys.gamma * 1e6, 0); % [uT]
rf_waveform2 = interp1(t_rf2, rf2, t_RRT2, 'linear', 'extrap');

%% Calculate a gradient waveform [mT/m]
grad_waveform1 = rf_waveform1 * 0;
grad_waveform2 = rf_waveform2 * 0;

%% Perform Bloch simulation
T1 = inf; % T1 relaxation time [sec]
T2 = inf; % T2 relaxation time [sec]

%--------------------------------------------------------------------------
% Calculate the range of off-resonance [Hz]
%--------------------------------------------------------------------------
Nf = 1000;
df_range = linspace(0, 1000, Nf); % [Hz]

%--------------------------------------------------------------------------
% Calculate the range of B1 [uT]
%--------------------------------------------------------------------------
N = 400;
% [G] * [T/1e4G] * [1e6uT/T] => *1e2 [uT]
B1_range = linspace(0, 1 * 1e2, N);

%--------------------------------------------------------------------------
% Calculate the initial magnetization
%--------------------------------------------------------------------------
mx0 = zeros(1, Nf, 'double');
my0 = zeros(1, Nf, 'double');
mz0 = ones(1, Nf, 'double');

%--------------------------------------------------------------------------
% Bloch simulation for a 4.7-ms BIR-4 pulse
%--------------------------------------------------------------------------
Mx1 = zeros(Nf, N, 'double');
My1 = zeros(Nf, N, 'double');
Mz1 = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    scale_factor = B1_range(idx) / b1_max;
    [mx,my,mz] = bloch(scale_factor * rf_waveform1 * 1e-2, grad_waveform1 * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx1(:,idx) = mx;
    My1(:,idx) = my;
    Mz1(:,idx) = mz;
end

mxy1 = Mx1 + 1j * My1;

%--------------------------------------------------------------------------
% Bloch simulation for a 5.7-ms BIR-4 pulse
%--------------------------------------------------------------------------
Mx2 = zeros(Nf, N, 'double');
My2 = zeros(Nf, N, 'double');
Mz2 = zeros(Nf, N, 'double');

for idx = 1:N
    start_tic = tic;
    fprintf('Performing Bloch simulation (%d/%d) ...', idx, N);
    % RF: [uT] * [T/1e6uT] * [1e4G/T] => *1e-2 [G]
    % GR: [mT/m] * [T/1e3mT] * [1e4G/T] * [m/1e2cm] => *1e-1 [G/cm]
    scale_factor = B1_range(idx) / b1_max;
    [mx,my,mz] = bloch(scale_factor * rf_waveform2 * 1e-2, grad_waveform2 * 1e-1, sys.rfRasterTime, T1, T2, df_range, 0, 0, mx0, my0, mz0);
    fprintf('done! (%6.4f sec)\n', toc(start_tic));
    Mx2(:,idx) = mx;
    My2(:,idx) = my;
    Mz2(:,idx) = mz;
end

mxy2 = Mx2 + 1j * My2;

%% Display Figure 7
[df_grid, B1_grid] = ndgrid(df_range, B1_range);

t1 = (-floor(N_RRT1/2):ceil(N_RRT1/2)-1).' * sys.rfRasterTime; % [sec]
t2 = (-floor(N_RRT2/2):ceil(N_RRT2/2)-1).' * sys.rfRasterTime; % [sec]

FontSize = 14;
figure('Color', 'w', 'Position', [-3 2 960 990]);
ax1 = subplot(3,2,[1,2]); 
hold on;
plot(t1 * 1e3, abs(rf_waveform1), 'LineWidth', 2, 'Color', [150 150 150] / 255);
plot(t2 * 1e3, abs(rf_waveform2), '--', 'LineWidth', 2, 'Color', 'k');
set(gca, 'XTickLabel', [], 'FontSize', FontSize, 'Box', 'on');
grid on;
xlim([-3 3]);
ylim([0 1.2]);
ylabel('$A(t)$ [unitless]', 'Interpreter', 'latex', 'FontSize', FontSize);
title('BIR-4', 'Interpreter', 'latex', 'FontSize', FontSize);

ax2 = subplot(3,2,[3,4]);
hold on;
plot(t1 * 1e3, angle(rf_waveform1), 'LineWidth', 2, 'Color', [150 150 150] / 255);
plot(t2 * 1e3, angle(rf_waveform2), '--', 'LineWidth', 2, 'Color', 'k');
set(gca, 'FontSize', FontSize, 'Box', 'on');
grid on;
xlim([-3 3]);
ylim([-4 4]);
xlabel('Time (ms)', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('$\phi(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', FontSize);

ax3 = subplot(3,2,5);
surf(df_grid, B1_grid * 1e-2, abs(mxy1), 'EdgeColor', 'none');
axis square;
view(90, -90);
colormap(turbo(256));
%hc = colorbar;
%ylabel(hc, '$|M_{xy}| (/M_0)$', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
xlabel('$\Delta f$ (Hz)', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel('$|B_1^+|$ (Gauss)', 'Interpreter', 'latex', 'FontSize', FontSize);
title(sprintf('BIR-4 %3.1f ms', Tp1 * 1e3), 'Interpreter', 'latex', 'FontSize', 16);

ax4 = subplot(3,2,6);
surf(df_grid, B1_grid * 1e-2, abs(mxy2), 'EdgeColor', 'none');
set(gca, 'XTickLabel', []);
axis square;
view(90, -90);
colormap(turbo(256));
hc = colorbar;
ylabel('$|B_1^+|$ (Gauss)', 'Interpreter', 'latex', 'FontSize', FontSize);
ylabel(hc, '$|M_{xy}| (/M_0)$', 'Interpreter', 'latex', 'FontSize', FontSize);
set(gca, 'FontSize', FontSize);
title(sprintf('BIR-4 %3.1f ms', Tp2 * 1e3), 'Interpreter', 'latex', 'FontSize', 16);

set(ax1, 'Position', [0.1300 0.7636 0.7750 0.1614]);
set(ax2, 'Position', [0.1300 0.5834 0.7750 0.1614]);

set(ax3, 'Position', [0.1300-0.01 0.1100+0.04 0.3347+0.1 0.2145+0.1]);
set(ax4, 'Position', [0.5703-0.06 0.1100+0.04 0.2458+0.1 0.2145+0.1]);

set(hc, 'Position', [0.9050-0.05 0.1100+0.04 0.0222 0.3145]);

export_fig(sprintf('Grissom_2014_JMR_figure7'), '-r300', '-tif');
