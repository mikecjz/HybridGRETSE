classdef TSE 
    
    properties
        
        dG=170e-6; %Gradient rise time
                
        tEx=0.4e-3; %Exitation RF time
        tRef=0.4e-3; %Refocusing RF Time
        
        fspR= 0.5; %Readout spoiler (crusher) area factor
        
        
        rfex_phase=pi/2; %Excitation RF phase
        rfref_phase=0; %Refocusing RF phase
        
        flipex=90*pi/180; %Exitation RF flip angle
        flipref=180*pi/180; %Refocusing RF flip angle
        
        
        %%%%%% The following properties are calculated in the class constructor %%%%%%%
        
        %Structure that defines the following properties for scan
        %parameters:
        %   fov : [xfov, yfov, zfov] in m
        %   Nx : Number of samples in the readout direction
        %   Ny : Number of samples in the phase encoding direction
        %   Nz : Number of samples in the partition direction
        %   nechos : TSE echo train length
        %   samplingTime : readout sampling time in s
        %   echoSpacing : Echo spacing time (or TE) in s
        scanParams = [];
        
        readoutTime = []; %Total readout gradient time: samplingTime + 2*adc_deadtime in s
        tExwd = []; % Excitation RF time with dead time in s
        tRefwd = []; % Refocusing RF time with dead time in s
        tSp = []; % Spoiler (crusher gradient) time in s
        tSpex = []; %Exitation spoiler (Used as Readout Prephase) time in s
        
        
        RFex = []; %RF excitation object (struct)
        RFref = []; %RF refocusing object (struct)
        
        
        GR3 = []; %Prephasing readout gradient
        GR5 = []; %Readout gradient spoiler (crusher) #1
        GR6 = []; %Readout gradient
        GR7 = []; %Readout fradient spoiler (crusher) #2
        
        ADC = []; %ADC object (struct)
        
        PE_area_steps = []; % Array of phase encoding (Ky) gradient areas m^-1
        Par_area_steps = []; % Array of partition encoding (Kz) gradient areas m^-1
        
    end
    
    methods
        
        %TSE object constructor
        function obj = TSE(system,scanParams)
            
            if nargin == 1 
                %In case scanParams are not specified, we use defaults
                newScanParams.fov = [256e-3, 256e-3, 256e-3]; % [xfov, yfov, zfov]
                newScanParams.Nx = 128;
                newScanParams.Ny = 128;
                newScanParams.Nz = 128;
                newScanParams.nechos = 60;
                newScanParams.samplingTime = 4e-3;
                newScanParams.echoSpacing = 6e-3;
                
                obj.scanParams = newScanParams;
                
            elseif nargin == 2
                obj.scanParams = scanParams;
            else
                error('TSE Constructor: Number of parameters is not correct. Usage: TSE(system,scanParams)')
                
            end
            
            %========== Calculate TSE object properties ========
            
            %Enforce RO and adc raster
            disp('Enforcing RO and adc raster, sampling Time might be slightly changed')
            obj.scanParams.samplingTime = ceil(obj.scanParams.samplingTime/obj.scanParams.Nx/system.adcRasterTime)...
                * system.adcRasterTime * obj.scanParams.Nx;% ADC raster
            
            % Gradient raster
            obj.scanParams.samplingTime = ceil(obj.scanParams.samplingTime./system.gradRasterTime) * system.gradRasterTime;
            
            obj.readoutTime = obj.scanParams.samplingTime + 2*system.adcDeadTime;
            
            obj.tExwd = obj.tEx + system.rfRingdownTime + system.rfDeadTime;
            obj.tRefwd = obj.tRef + system.rfRingdownTime + system.rfDeadTime;
            obj.tSp = 0.5 * (obj.scanParams.echoSpacing - obj.readoutTime - obj.tRefwd);
            obj.tSpex = 0.5 * (obj.scanParams.echoSpacing - obj.tExwd - obj.tRefwd);
            
            
            [obj.RFex, delayEx] = mr.makeBlockPulse(obj.flipex,system,'PhaseOffset',obj.rfex_phase,'Duration',obj.tEx);
            [obj.RFref, delayRef] = mr.makeBlockPulse(obj.flipref,system,'PhaseOffset',obj.rfref_phase,'Duration',obj.tRef);
            
            
            
            fov = obj.scanParams.fov;
            Nx = obj.scanParams.Nx;
            Ny = obj.scanParams.Ny;
            Nz = obj.scanParams.Nz;
            
            deltaKx=1/fov(1);
            kxWidth = Nx*deltaKx;
            
            deltaKy = 1/fov(2);
            deltaKz = 1/fov(3);
            pe_steps = (1:(Ny))-0.5*Ny-1;
            par_steps = (1:(Nz))-0.5*Nz-1;
            
            obj.PE_area_steps = pe_steps*deltaKy;
            obj.Par_area_steps = par_steps*deltaKz;
            
            
            %Preliminary Readout acuisition gradient and spoiler gradient.
            %Used for later calculation of actual gradients
            GRacq = mr.makeTrapezoid('x',system,'FlatArea',kxWidth,'FlatTime',obj.readoutTime,'riseTime',obj.dG);
            GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area*obj.fspR,'duration',obj.tSp,'riseTime',obj.dG);
            
            AGRspr=GRspr.area;%GRacq.area/2*fspR;
            AGRpreph = GRacq.area/2+AGRspr;%GRacq.area*(1+fspR)/2;
            GRpreph = mr.makeTrapezoid('x',system,'Area',AGRpreph,'duration',obj.tSpex,'riseTime',obj.dG); %Prephase readout gradient
            
            
            %Readout objects calculation
            obj.GR3 = GRpreph;
            
            GR5times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
            GR5amp=[0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
            obj.GR5 = mr.makeExtendedTrapezoid('x','times',GR5times,'amplitudes',GR5amp);
            
            GR6times=[0 obj.readoutTime];
            GR6amp=[GRacq.amplitude GRacq.amplitude];
            obj.GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);
            
            GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
            GR7amp=[GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
            obj.GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);
            
            obj.ADC = mr.makeAdc(obj.scanParams.Nx,'Duration',obj.scanParams.samplingTime, 'Delay', system.adcDeadTime);
            obj.ADC.phaseOffset = pi/2; %Add ADC offset to match phase with GRE
            %Readout objects calculation END
            
            
            
            %========== Calculate TSE object properties END ========

        end
        
        %Append a TSE block to the seq object
        % INPUTS:
        %  -(obj: reference to this TSE object itself)
        %  - seq : The sequence object to be appended to
        %  - pe_indices: [ky, kz] Array of kspace locations. size(pe_indices,1) should be equal to nechos
        %  - system: system object (struct) keeping track of heardware limitations
        % OUTPUTS:
        %  - seq: appended seq object
        
        function seq = AppendTSEBlock(obj,seq,system,pe_indices)
           
            nechos = obj.scanParams.nechos;
            
            if size(pe_indices,1) ~= nechos
                error(['TSE Block append: size of pe_indices is not consistant with number of echos ',...
                    '(',num2str(nechos),')'])
            end
            
            ky_array = pe_indices(:,1);
            kz_array = pe_indices(:,2);
            
            PE_area_array = obj.PE_area_steps(ky_array);
            Par_area_array = obj.Par_area_steps(kz_array);
            
            %A spoil gradient in PE Gradient
            tPSpoil = 3000e-6;
            GP_spoil = mr.makeTrapezoid('y',system,'duration',tPSpoil,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
            seq.addBlock(GP_spoil);
            
            %Add Excitation RF pulse and prephase readout gradient
            seq.addBlock(obj.RFex);
            seq.addBlock(obj.GR3);
            
            %Add refocusing TSE train
            for iEcho = 1:nechos
                
                %Calculate phase and partition gradiatent from specified k-space location
                GPpre = mr.makeTrapezoid('y',system,'Area',PE_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                GPrew = mr.makeTrapezoid('y',system,'Area',-PE_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                
                GSpre = mr.makeTrapezoid('z',system,'Area',Par_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                GSrew = mr.makeTrapezoid('z',system,'Area',-Par_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                
                %Add refocusing RF pulse and subsequent gradients
                seq.addBlock(obj.RFref);
                
                seq.addBlock(obj.GR5,GPpre,GSpre);
                seq.addBlock(obj.GR6,obj.ADC);
                
                seq.addBlock(obj.GR7,GPrew,GSrew);
                
            end
            
            
        end
        
        %Only FID, No Gradients
        function seq = AppendTSEFIDBlock(obj,seq,system,pe_indices)
           
            nechos = obj.scanParams.nechos;
            
            if size(pe_indices,1) ~= nechos
                error(['TSE Block append: size of pe_indices is not consistant with number of echos ',...
                    '(',num2str(nechos),')'])
            end
            
            ky_array = pe_indices(:,1);
            kz_array = pe_indices(:,2);
            
            PE_area_array = obj.PE_area_steps(ky_array);
            Par_area_array = obj.Par_area_steps(kz_array);
            
            %A spoil gradient in PE Gradient
            tPSpoil = 3000e-6;
            GP_spoil = mr.makeTrapezoid('y',system,'duration',tPSpoil,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
%             seq.addBlock(GP_spoil);

            GR3Delay = mr.makeDelay(mr.calcDuration(obj.GR3));
            GR5Delay = mr.makeDelay(mr.calcDuration(obj.GR5));
            GR6Delay = mr.makeDelay(mr.calcDuration(obj.GR6));
            GR7Delay = mr.makeDelay(mr.calcDuration(obj.GR7));
            
            %Add Excitation RF pulse and prephase readout gradient
            seq.addBlock(obj.RFex);
            seq.addBlock(GR3Delay);
            
            
            
            %Add refocusing TSE train
            for iEcho = 1:nechos
                
                
                
                %Add refocusing RF pulse and subsequent gradients
                seq.addBlock(obj.RFref);
                
                seq.addBlock(GR5Delay);
                seq.addBlock(GR6Delay,obj.ADC);
                
                seq.addBlock(GR7Delay);
                
            end
            
            
        end
        
        %Append a VFL TSE block to the seq object
        % INPUTS:
        %  -(obj: reference to this TSE object itself)
        %  - seq : The sequence object to be appended to
        %  - pe_indices: [ky, kz] Array of kspace locations. size(pe_indices,1) should be equal to nechos
        %  - system: system object (struct) keeping track of heardware limitations
        % OUTPUTS:
        %  - seq: appended seq object
        
        function seq = AppendTSE_VFL_Block(obj,seq,system,pe_indices)
           
            nechos = obj.scanParams.nechos;
            
            if size(pe_indices,1) ~= nechos
                error(['TSE Block append: size of pe_indices is not consistant with number of echos ',...
                    '(',num2str(nechos),')'])
            end
            
            ky_array = pe_indices(:,1);
            kz_array = pe_indices(:,2);
            
            PE_area_array = obj.PE_area_steps(ky_array);
            Par_area_array = obj.Par_area_steps(kz_array);
            
            faArray = TSE.CalculateVFLAngles(nechos,'generic');
            
            %A spoil gradient in PE Gradient
            tPSpoil = 3000e-6;
            GP_spoil = mr.makeTrapezoid('y',system,'duration',tPSpoil,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
            seq.addBlock(GP_spoil);
            
            %Add Excitation RF pulse and prephase readout gradient
            seq.addBlock(obj.RFex);
            seq.addBlock(obj.GR3);
           
            
            %Add refocusing TSE train
            for iEcho = 1:nechos
                
                tmpRFref = mr.makeBlockPulse(faArray(iEcho),system,'PhaseOffset',obj.rfref_phase,'Duration',obj.tRef);
                
                %Calculate phase and partition gradiatent from specified k-space location
                GPpre = mr.makeTrapezoid('y',system,'Area',PE_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                GPrew = mr.makeTrapezoid('y',system,'Area',-PE_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                
                GSpre = mr.makeTrapezoid('z',system,'Area',Par_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                GSrew = mr.makeTrapezoid('z',system,'Area',-Par_area_array(iEcho),'Duration',obj.tSp,'riseTime',obj.dG);
                
                %Add refocusing RF pulse and subsequent gradients
                seq.addBlock(tmpRFref);
                
                seq.addBlock(obj.GR5,GPpre,GSpre);
                seq.addBlock(obj.GR6,obj.ADC);
                
                seq.addBlock(obj.GR7,GPrew,GSrew);
                
            end
            
            
        end
        
        
        
    end
    
    methods(Static)
        %Outputs flip angles array in radian
        function faArray = CalculateVFLAngles(nEchos, mode)
            
            if strcmp(mode,'generic')
                nDecrease = 15;
                VFLBounds = [120,150];% Degree
                VFLRange = abs(VFLBounds(1)-VFLBounds(2));
                
                tempX = linspace(0,5,nDecrease);
                
                faDecrease = round(exp(-tempX)*VFLRange+min(VFLBounds(:)));
                
                faArray = min(VFLBounds(:))*ones(1,nEchos);
                faArray(1:nDecrease) = faDecrease;
                
                faArray = faArray * pi/180;
                
            end
            
        end
    end
    
end