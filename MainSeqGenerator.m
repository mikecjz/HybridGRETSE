%% Include pulseq in path

addpath('pulseq/matlab')

%%

system = mr.opts('MaxGrad', 35/sqrt(3),... 
                'GradUnit', 'mT/m', ...
                'MaxSlew', 190/sqrt(3),...
                'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 100e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6,...
                'B0', 0.55);
            
%% Global parameters
fov = [256e-3, 256e-3, 256e-3];
Nx = 128;
Ny = 128;
Nz = 128;

%% TSE parameters
TSE_scanParams.fov = fov;
TSE_scanParams.Nx = Nx;
TSE_scanParams.Ny = Ny;
TSE_scanParams.Nz = Nz;
TSE_scanParams.nechos = 60;
TSE_scanParams.samplingTime = 6.4e-3;
TSE_scanParams.echoSpacing  = 12e-3;

%% GRE parameters
GRE_scanParams.fov = fov;
GRE_scanParams.Nx = Nx;
GRE_scanParams.Ny = Ny;
GRE_scanParams.Nz = Nz;
GRE_scanParams.nechos = [150,150];
GRE_scanParams.samplingTime = 3e-3;
GRE_scanParams.echoSpacing = 6e-3;

%% Caculate Spiral In trajectory
nav_interval = 10;
sp_k = 10;%spiral parameterk

% [ky,kz] = ROCK(1:100000,nav_interval,sp_k,[Nz,Ny],0.5e-2);
load('traj.mat')
ky = ky(:);
kz = kz(:);


ky(:) = 65*ones(100000,1);
kz(:) = 65*ones(100000,1);
%% Define sequence

clear seq

% Define TR Gap in needed
TRfill=10e-3;
% round to gradient raster
TRfill=system.gradRasterTime * round(TRfill / system.gradRasterTime);
if TRfill<0, TRfill=1e-3; 
    disp(strcat('Warning!!! TR too short, adapted to include all slices to : ',num2str(1000*Nslices*(tETrain+TRfill)),' ms')); 
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms')); 
end
delayTR = mr.makeDelay(TRfill);

seq=mr.Sequence(system);

MT_TSE = TSE(system,TSE_scanParams);
MT_GRE = GRE(system,GRE_scanParams);

tseEchos = TSE_scanParams.nechos;
greEchos = sum(GRE_scanParams.nechos(:));
totalEchos = tseEchos + greEchos;

for iTR = 1:10
    disp(iTR)
    gre_indices = (1:greEchos) + (iTR-1) * totalEchos;
    pe_indices_gre = [ky(gre_indices),kz(gre_indices)];
    
    tse_indices = ((greEchos+1):(greEchos+tseEchos)) + (iTR-1) * totalEchos;
    pe_indices_tse = [ky(tse_indices),kz(tse_indices)];
    
    seq = MT_GRE.AppendGREBlock(seq,system,pe_indices_gre);
    seq.addBlock(delayTR);
    seq = MT_TSE.AppendTSEBlock(seq,system,pe_indices_tse);
    
   
    
end

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
%% Write Sequence
save_dir = '/mnt/radnas1/Junzhou/Scanner_Bins/PulseSeq';

seq.write(fullfile(save_dir,'GRETSE_nav.seq'))
% save('traj.mat','ky','kz')