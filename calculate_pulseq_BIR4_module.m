function rf_bir4 = calculate_pulseq_BIR4_module(T_seg, b1_max, dw_max, zeta, kappa, beta, sys)
% calculate_pulseq_BIR4_module -- Calculate a BIR-4 module
%
%  Usage
%    rf_bir4 = calculate_pulseq_BIR4_module(T_seg, b1_max, dw_max, zeta, kappa, beta, sys)
%  Inputs
%    T_seg       duration of one pulse segment [sec]
%    b1_max      maximum RF amplitude [uT]
%    dw_max      maximum frequency sweep [Hz]
%    zeta        constant in the amplitude function [rad]
%    kappa       constant in the frequency/phase function [rad]
%    beta        flip angle [degree]
%    sys         Pulseq "system"
%  Outputs
%    rf_bir4     Pulseq RF event
%
%  Description
%    BIR-4 based B1 and B0 insensitive velocity selective pulse train described
%    by Wong and Guo 2010 ISMRM
%
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/08/2022, Last modified: 09/08/2022

%% BIR-4 module
%--------------------------------------------------------------------------
%
%          RF       RF       RF       RF 
%       segment1 segment2 segment3 segment4
%      |<------>|<------>|<------>|<------>|
%       _______   _______|_______   _______         
%      |       \ /       |       \ /       |
%      |        |        |        |        |
%      |        |        |        |        |
% o----+--------+--------+--------+--------+-> t
% |<-->|<------>|<------>|<------>|<------>|
% deadTime T_seg| T_seg    T_seg  | T_seg
%--------------------------------------------------------------------------

%% Recalculate T_seg and Tp
%--------------------------------------------------------------------------
% Set the number of samples in one segment as an even number
%--------------------------------------------------------------------------
N_seg = floor(T_seg / sys.rfRasterTime);
if mod(N_seg,2) == 1 % odd
    N_seg = N_seg + 1;
end
T_seg = N_seg * sys.rfRasterTime;

%--------------------------------------------------------------------------
% Calculate the duration of all segments [sec]
%--------------------------------------------------------------------------
Tp = 4 * T_seg;

%% Calculate a BIR-4 module
rf_samples = N_seg * 4; % number of RF samples
t = ((0:rf_samples-1).' + 0.5) * sys.rfRasterTime; % RRT: RF RASTER TIME

%--------------------------------------------------------------------------
% Calculate the maximum RF amplitude in [Hz]
%--------------------------------------------------------------------------
% [uT] * [Hz/T] * [T/1e6uT] => * 1e-6 [Hz]
w1_max = b1_max * sys.gamma * 1e-6; % [Hz]

%--------------------------------------------------------------------------
% Define dphi1 and dphi2
%--------------------------------------------------------------------------
dphi1 =  (180 + beta / 2) * pi / 180; % [rad]
dphi2 = -(180 + beta / 2) * pi / 180; % [rad]

%--------------------------------------------------------------------------
% Calculate phi_max
%--------------------------------------------------------------------------
tan_kappa = tan(kappa);
% [Hz] * [2pi rad/cycle] * [sec] / ([rad] * [unitless]) => [rad]
phi_max = -(dw_max * 2 * pi) * Tp / (kappa * tan_kappa) * log(cos(kappa)); % [rad]

%--------------------------------------------------------------------------
% Define function handles for amplitude and phase functions
%--------------------------------------------------------------------------
w1  = @(t) w1_max * tanh(zeta * (1 - 4 * t / Tp)); % [Hz]
phi = @(t) phi_max / 4 - (dw_max * 2 * pi) * Tp / (4 * kappa * tan_kappa) * log(cos(4 * kappa * t / Tp) / cos(kappa)); % [rad]

%--------------------------------------------------------------------------
% Segment 1 (0 < t <= 0.25 * Tp)
%--------------------------------------------------------------------------
segment1_range = (t <= 0.25 * Tp);
t1 = t(segment1_range);
am_segment1 = w1(t1);
pm_segment1 = phi(t1);

%--------------------------------------------------------------------------
% Segment 2 (0.25 * Tp < t <= 0.5 * Tp)
%--------------------------------------------------------------------------
segment2_range = (t > 0.25 * Tp) & (t <= 0.5 * Tp);
t2 = t(segment2_range);
am_segment2 = w1(0.5 * Tp - t2);
pm_segment2 = phi(0.5 * Tp - t2) + dphi1;

%--------------------------------------------------------------------------
% Segment 3 (0.5 * Tp < t <= 0.75 * Tp)
%--------------------------------------------------------------------------
segment3_range = (t > 0.5 * Tp) & (t <= 0.75 * Tp);
t3 = t(segment3_range);
am_segment3 = w1(t3 - 0.5 * Tp);
pm_segment3 = phi(t3 - 0.5 * Tp) + dphi1;

%--------------------------------------------------------------------------
% Segment 4 (0.75 * Tp < t <= Tp)
%--------------------------------------------------------------------------
segment4_range = (t > 0.75 * Tp) & (t <= Tp);
t4 = t(segment4_range);
am_segment4 = w1(Tp - t4);
pm_segment4 = phi(Tp - t4) + dphi1 + dphi2;

%--------------------------------------------------------------------------
% Combine all segments to form BIR-4
%--------------------------------------------------------------------------
am = cat(1, am_segment1, am_segment2, am_segment3, am_segment4);
pm = cat(1, pm_segment1, pm_segment2, pm_segment3, pm_segment4);
rf_shape = am .* exp(1j * pm);

%% Create an RF event
rf_bir4 = mr.makeArbitraryRf(rf_shape * 0, pi, 'system', sys);
rf_bir4.signal = rf_shape;

end