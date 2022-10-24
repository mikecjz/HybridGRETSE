system = mr.opts('MaxGrad', 35,...
    'GradUnit', 'mT/m', ...
    'MaxSlew', 190/sqrt(3),...
    'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6,...
    'B0', 0.55);

%% Define parameters for a BIR-4 VS preparation
zeta   = 10;                       % constant in the amplitude function [rad]
beta   = 45;                       % flip angle [degree]
b1_max = 10;                        % maximum RF amplitude [uT]

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
rf1_bir4 = calculate_pulseq_BIR4_module(T_seg1, b1_max, dw_max1, zeta, kappa1, beta, system);

tRFref = mr.calcDuration(rf1_bir4) + 16e-6;

rf_ex  = mr.makeBlockPulse(90*pi/180,'system', system,'PhaseOffset',90*pi/180,'Duration',rf1_bir4.shape_dur+ 16e-6);





clear seq

seq=mr.Sequence(system);

samplingTime = 6.4e-3 + 2*system.adcDeadTime;

ADC = mr.makeAdc(128,'Duration',6.4e-3 , 'Delay', system.adcDeadTime);

ADCDelay = mr.makeDelay(3.6e-3);

ES = samplingTime + 2 * ADCDelay.delay + tRFref;

EXDelay = mr.makeDelay( 0.5*ES - tRFref);

delayTR = mr.makeDelay(4.0);

for iTR = 1:5
    seq.addBlock(rf_ex);
    seq.addBlock(EXDelay);
   for iEcho = 1:60
       seq.addBlock(rf1_bir4,mr.makeDelay(tRFref));
      
       seq.addBlock(ADCDelay);
       seq.addBlock(mr.makeDelay(samplingTime),ADC);
       seq.addBlock(ADCDelay);
   end
   seq.addBlock(delayTR);
end


[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
%%
save_dir = '/mnt/radnas1/Junzhou/Scanner_Bins/PulseSeq';

seq.write(fullfile(save_dir,'TSEFID12ES_BIR4.seq'))
% save('traj.mat','ky','kz')