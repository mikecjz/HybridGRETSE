%% Create a TSE sequence and export for execution
% 
% The |Sequence| class provides functionality to create magnetic
% resonance sequences (MRI or NMR) from basic building blocks.
%
% This provides an implementation of the open file format for MR sequences
% described here: http://pulseq.github.io/specification.pdf
%
% This example performs the following steps:
% 
% # Create slice selective RF pulse for imaging.
% # Create readout gradient and phase encode strategy.
% # Loop through phase encoding and generate sequence blocks.
% # Write the sequence to an open file format suitable for execution on a
% scanner.
% 
%   Juergen Hennig <juergen.hennig@uniklinik-freiburg.de>
%   Maxim Zaitsev  <maxim.zaitsev@uniklinik-freiburg.de>
 

%% Instantiation and gradient limits
% The system gradient limits can be specified in various units _mT/m_,
% _Hz/cm_, or _Hz/m_. However the limits will be stored internally in units
% of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecificied hardware
% parameters will be assigned default values.

dG=250e-6;
system = mr.opts('MaxGrad', 35/sqrt(3),... 
                'GradUnit', 'mT/m', ...
                'MaxSlew', 190/sqrt(3),...
                'SlewUnit', 'T/m/s', ...
                'rfRingdownTime', 100e-6, ...
                'rfDeadTime', 100e-6, ...
                'adcDeadTime', 10e-6,...
                'B0', 0.55);

%%
% A new sequence object is created by calling the class constructor.
seq=mr.Sequence(system);


%% Sequence events
% Some sequence parameters are defined using standard MATLAB variables
fov=256e-3;
Nx=128; Ny=128; necho=60; Nslices=1;
rflip=180;
if (numel(rflip)==1), rflip=rflip+zeros([1 necho]); end
sliceThickness=5e-3;
TE=12e-3; TR=2000e-3;
TEeff=60e-3;
k0=round(TEeff/TE);
PEtype='linear';

samplingTime= 6.4e-3;
%Enforce GO raster
samplingTime = round(samplingTime/Nx/system.gradRasterTime) * system.gradRasterTime * Nx;

readoutTime = samplingTime + 2*system.adcDeadTime;
tEx=2.5e-3; 
tExwd=tEx+system.rfRingdownTime+system.rfDeadTime;
tRef=2.0e-3; 
tRefwd=tRef+system.rfRingdownTime+system.rfDeadTime;
tSp=0.5*(TE-readoutTime-tRefwd);
tSpex=0.5*(TE-tExwd-tRefwd);
fspR=1.0;
fspS=0.5;

rfex_phase=pi/2; % MZ: we need to maintain these as variables because we will overwrtite phase offsets for multiple slice positions
rfref_phase=0;

%%
%%% Base gradients
%%% Slice selection
% Key concepts in the sequence description are *blocks* and *events*.
% Blocks describe a group of events that are executed simultaneously. This
% hierarchical structure means that one event can be used in multiple
% blocks, a common occurrence in MR sequences, particularly in imaging
% sequences. 
%
% First, the slice selective RF pulses (and corresponding slice gradient)
% are generated using the |makeSincPulse| function.
% Gradients are recalculated such that their flattime covers the pulse plus
% the rfdead- and rfringdown- times.
%
flipex=90*pi/180;
[rfex,delayEx] = mr.makeBlockPulse(flipex,system,'PhaseOffset',rfex_phase,'Duration',tEx);


flipref=rflip(1)*pi/180;
[rfref,delayRef] = mr.makeBlockPulse(flipref,system,'PhaseOffset',rfref_phase,'Duration',tRef);

%%
%%% Readout gradient
% To define the remaining encoding gradients we need to calculate the
% $k$-space sampling. The Fourier relationship
%
% $$\Delta k = \frac{1}{FOV}$$
% 
% Therefore the area of the readout gradient is $n\Delta k$.
deltak=1/fov;
kWidth = Nx*deltak;

GRacq = mr.makeTrapezoid('x',system,'FlatArea',kWidth,'FlatTime',readoutTime,'riseTime',dG);
adc = mr.makeAdc(Nx,'Duration',samplingTime, 'Delay', system.adcDeadTime);%,'Delay',GRacq.riseTime);
GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area*fspR,'duration',tSp,'riseTime',dG);
GRspex = mr.makeTrapezoid('x',system,'area',GRacq.area*(1+fspR),'duration',tSpex,'riseTime',dG);


AGRspr=GRspr.area;%GRacq.area/2*fspR;
AGRpreph = GRacq.area/2+AGRspr;%GRacq.area*(1+fspR)/2;
GRpreph = mr.makeTrapezoid('x',system,'Area',AGRpreph,'duration',tSpex,'riseTime',dG);



%%
%%% Phase encoding
% To move the $k$-space trajectory away from 0 prior to the readout a
% prephasing gradient must be used. Furthermore rephasing of the slice
% select gradient is required.

%[PEorder,Ny] = myTSE_PEorder(Ny,necho,k0,PEtype);
nex=floor(Ny/necho);
pe_steps=(1:(Ny))-0.5*Ny-1;
phaseAreas = pe_steps*deltak;


nav_interval = 10;
sp_k = 10;%spiral parameterk

[ky,kz] = ROCK(1:16800,nav_interval,sp_k,[128,128],0.5e-2);

ky(:) = 65;
kz(:) = 65;

%% split gradients and recombine into blocks


% and now the readout gradient....

GR3=GRpreph;%GRspex;

GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);

GR6times=[0 readoutTime];
GR6amp=[GRacq.amplitude GRacq.amplitude];
GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);

GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR7amp=[GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);


% and filltimes



TRfill=1.8;
% round to gradient raster
TRfill=system.gradRasterTime * round(TRfill / system.gradRasterTime);
if TRfill<0, TRfill=1e-3; 
    disp(strcat('Warning!!! TR too short, adapted to include all slices to : ',num2str(1000*Nslices*(tETrain+TRfill)),' ms')); 
else
    disp(strcat('TRfill : ',num2str(1000*TRfill),' ms')); 
end
delayTR = mr.makeDelay(TRfill);

%% Define sequence blocks
% Next, the blocks are put together to form the sequence
iLinesIndex = 0;
for iEchoTrains = 1:10
    seq.addBlock(rfex);
    seq.addBlock(GR3);
    for iEcho = 1:necho
        iLinesIndex = iLinesIndex+1;
        GPpre = mr.makeTrapezoid('y',system,'Area',phaseAreas(ky(iLinesIndex)),'Duration',tSp,'riseTime',dG);
        GPrew = mr.makeTrapezoid('y',system,'Area',-phaseAreas(ky(iLinesIndex)),'Duration',tSp,'riseTime',dG);
        
        GSpre = mr.makeTrapezoid('z',system,'Area',phaseAreas(kz(iLinesIndex)),'Duration',tSp,'riseTime',dG);
        GSrew = mr.makeTrapezoid('z',system,'Area',-phaseAreas(kz(iLinesIndex)),'Duration',tSp,'riseTime',dG);
        
        seq.addBlock(rfref);
        seq.addBlock(GR5,GPpre,GSpre);
        seq.addBlock(GR6,adc);
        
        seq.addBlock(GR7,GPrew,GSrew);
    end
    seq.addBlock(delayTR);
    
end

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
% 
% % plot k-spaces
% figure; plot(t_ktraj,ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
% figure; plot(ktraj(1,:),ktraj(2,:),'b',...
%              ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% title('2D k-space');


%% Write to file

% The sequence is written to file in compressed form according to the file
% format specification using the |write| method.
seq.write('tse_spin.seq')

%%
% Display the first few lines of the output file
% s=fileread('myTSE.seq');
% disp(s(1:300))
% seq.plot();

% seq.install('siemens');

