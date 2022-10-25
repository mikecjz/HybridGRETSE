classdef GRE 
    
    properties
        
        dG=170e-6; %Gradient rise time
                
        tEx=0.4e-3; %Exitation RF time
        
        fspR=0.5; %Readout spoiler (crusher) area factor
        
      
        
        flipex=[5,1] * pi/180; %Exitation RF flip angle
        
        
        
        %%%%%% The following properties are calculated in the class constructor %%%%%%%
        
        %Structure that defines the following properties for scan
        %parameters:
        %   fov : [xfov, yfov, zfov] in m
        %   Nx : Number of samples in the readout direction
        %   Ny : Number of samples in the phase encoding direction
        %   Nz : Number of samples in the partition direction
        %   nechos : GRE echo train length of each segments e.g. [150,150]
        %   samplingTime : readout sampling time in s
        %   echoSpacing : Echo spacing time (or TE) in s
        scanParams = [];
        
        readoutTime = []; %Total readout gradient time: samplingTime + 2*adc_deadtime in s
        tExwd = []; % Excitation RF time with dead time in s
        
        RFex1 = []; %First Set of RF object (default 5 degree)
        RFex2 = []; %Second Set of RF object (default 1 degree)
        
        tSp = []; % Spoiler (crusher gradient) time in s
        tPreph = []; % Prephasing gradient time in s
        
        
        GR3 = []; %Prephasing readout gradient
        
        GR6 = []; %Readout gradient
        GR7 = []; %Readout gradient spoiler 
        
        ADC = []; %ADC object (struct)
        
        PE_area_steps = []; % Array of phase encoding (Ky) gradient areas m^-1
        Par_area_steps = []; % Array of partition encoding (Kz) gradient areas m^-1
        
    end
    
    methods
        
        %GRE object constructor
        function obj = GRE(system,scanParams)
            
            if nargin == 1 
                %In case scanParams are not specified, we use defaults
                newScanParams.fov = [256e-3, 256e-3, 256e-3]; % [xfov, yfov, zfov]
                newScanParams.Nx = 128;
                newScanParams.Ny = 128;
                newScanParams.Nz = 128;
                newScanParams.nechos = [150,150];
                
                newScanParams.samplingTime = 4e-3;
                newScanParams.echoSpacing = 6e-3;
                
                obj.scanParams = newScanParams;
                
            elseif nargin == 2
                obj.scanParams = scanParams;
            else
                error('GRE Constructor: Number of parameters is not correct. Usage: GRE(system,scanParams)')
                
            end
            
            %========== Calculate TSE object properties ========
            
            %Enforce GO raster
            disp('Enforcing RO raster, sampling Time might be slightly changed')
            obj.scanParams.samplingTime = round(obj.scanParams.samplingTime/obj.scanParams.Nx/system.gradRasterTime)...
                * system.gradRasterTime * obj.scanParams.Nx;
            
            obj.readoutTime = obj.scanParams.samplingTime + 2*system.adcDeadTime;
            
            obj.tExwd = obj.tEx + system.rfRingdownTime + system.rfDeadTime;
            
            obj.RFex1 = mr.makeBlockPulse(obj.flipex(1),system,'PhaseOffset',0,'Duration',obj.tEx);
            obj.RFex2 = mr.makeBlockPulse(obj.flipex(2),system,'PhaseOffset',0,'Duration',obj.tEx);


            
            
            
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
            GRacq = mr.makeTrapezoid('x',system,'FlatArea',kxWidth,'FlatTime',obj.readoutTime);%,'riseTime',obj.dG);
            ROGradientTime = obj.readoutTime + GRacq.riseTime + GRacq.fallTime;
            
            obj.tSp = 2/4 * (obj.scanParams.echoSpacing - ROGradientTime - obj.tExwd);
            obj.tPreph = 2/4 * (obj.scanParams.echoSpacing - ROGradientTime - obj.tExwd);
            
            
            GRspr = mr.makeTrapezoid('x',system,'area',GRacq.area,'duration',obj.tSp);%,'riseTime',obj.dG);
            
            AGRspr=GRspr.area; %GRacq.area/2*fspR;
            
            AGRpreph = -1* GRacq.area/2;%GRacq.area*(1+fspR)/2;
            
            GRpreph = mr.makeTrapezoid('x',system,'Area',AGRpreph,'duration',obj.tPreph);%,'riseTime',obj.dG); %Prephase readout gradient
            
            
            %Readout objects calculation
            obj.GR3 = GRpreph;

            
            GR6times=[0, GRacq.riseTime, GRacq.riseTime + obj.readoutTime];
            GR6amp=[0, GRacq.amplitude, GRacq.amplitude];
            obj.GR6 = mr.makeExtendedTrapezoid('x','times',GR6times,'amplitudes',GR6amp);
            
            GR7times=[0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
            GR7amp=[GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
            obj.GR7 = mr.makeExtendedTrapezoid('x','times',GR7times,'amplitudes',GR7amp);
            
            obj.ADC = mr.makeAdc(obj.scanParams.Nx,'Duration',obj.scanParams.samplingTime, 'Delay', system.adcDeadTime + GRacq.riseTime);
            %Readout objects calculation END
            
            
            
            %========== Calculate TSE object properties END ========

        end
        
        %Append a GRE block to the seq object
        % INPUTS:
        %  -(obj: reference to this TSE object itself)
        %  - seq : The sequence object to be appended to
        %  - system: system object (struct) keeping track of heardware limitations
        %  - pe_indices: [ky, kz] Array of kspace locations. size(pe_indices,1) should be equal to nechos
        % OUTPUTS:
        %  - seq: appended seq object
        
        function seq = AppendGREBlock(obj,seq,system,pe_indices)
           
            totalEchos = sum(obj.scanParams.nechos(:));
            
            if size(pe_indices,1) ~= totalEchos
                error(['GRE Block append: size of pe_indices is not consistant with number of echos ',...
                    '(',num2str(totalEchos),')'])
            end
            
            ky_array = pe_indices(:,1);
            kz_array = pe_indices(:,2);
            
            PE_area_array = obj.PE_area_steps(ky_array);
            Par_area_array = obj.Par_area_steps(kz_array);
            
            
            
            % Add SR
            seq = obj.AppendSRBlock(seq,system);
            
            %Add GRE train
            for iEcho = 1:totalEchos
                
                % Define current RF for the 1st or 2nd GRE segment
                if iEcho <= obj.scanParams.nechos(1)
                    GRErf = obj.RFex1;
                else
                    GRErf = obj.RFex2;
                end
                
                %Adjust RF and ADC Phase (RF spoiling)
                GRErf.phaseOffset = mod(117*(iEcho^2+iEcho+2)*pi/180,2*pi);
                obj.ADC.phaseOffset = GRErf.phaseOffset;
                
                %Add Excitation RF pulse 
                seq.addBlock(GRErf);
                
                
                %Calculate phase and partition gradiatent from specified k-space location
                GPpre = mr.makeTrapezoid('y',system,'Area',PE_area_array(iEcho),'Duration',obj.tPreph);%,'riseTime',obj.dG);
                GPrew = mr.makeTrapezoid('y',system,'Area',-PE_area_array(iEcho),'Duration',obj.tPreph);%,'riseTime',obj.dG);
                
                GSpre = mr.makeTrapezoid('z',system,'Area',Par_area_array(iEcho),'Duration',obj.tPreph);%,'riseTime',obj.dG);
                GSrew = mr.makeTrapezoid('z',system,'Area',-Par_area_array(iEcho),'Duration',obj.tPreph);%,'riseTime',obj.dG);
                
                seq.addBlock(obj.GR3, GPpre, GSpre); %Pre phase Gradients
                
                seq.addBlock(obj.GR6,obj.ADC);
                
                seq.addBlock(obj.GR7,GPrew,GSrew);
                
            end
            
            
        end
        
        %Append a Saturation Recovery block to the seq object
        % INPUTS:
        %  -(obj: reference to this TSE object itself)
        %  - seq : The sequence object to be appended to- pe_indices: [ky, kz] Array of kspace locations. size(pe_indices,1) should be equal to nechos
        %  - system: system object (struct) keeping track of heardware limitations
        % OUTPUTS:
        %  - seq: appended seq object
        function seq = AppendSRBlock(obj,seq,system)
            % SR has 3 90 RF pulses. This makes 4 segments
            %      90     90     90
            %       |      |      |     
            %       |      |      |     
            % SEG 1 | SEG2 | SEG3 | SEG4
            
            % Define gradient segments length
            
            tSeg1 = 1000e-6;
            tSeg2 = 8800e-6;
            tSeg3 = 5800e-6;
            tSeg4 = 2550e-6;
            
            % Define RF
            RF90 = mr.makeBlockPulse(pi/2,system,'PhaseOffset',0,'Duration',obj.tEx);
            
           
            %Seg1 Gradients
            GR_seg1 = mr.makeTrapezoid('x',system,'duration',tSeg1,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
            GS_seg1 = mr.makeTrapezoid('z',system,'duration',tSeg1,'amplitude',-1* system.maxGrad);%,'riseTime',obj.dG);
            
            %Seg2 Gradients
            GR_seg2 = mr.makeTrapezoid('x',system,'duration',tSeg2,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
            GS_seg2 = mr.makeTrapezoid('z',system,'duration',tSeg2,'amplitude',system.maxGrad);%,'riseTime',obj.dG);
            
            %Seg3 Gradients
            GR_seg3 = mr.makeTrapezoid('x',system,'duration',tSeg3,'amplitude',-1*system.maxGrad);%,'riseTime',obj.dG);
            GS_seg3 = mr.makeTrapezoid('z',system,'duration',tSeg3,'amplitude',-1*system.maxGrad);%,'riseTime',obj.dG);
            
            %Seg4 Gradients
            GP_seg4 = mr.makeTrapezoid('y',system,'duration',tSeg4,'amplitude', system.maxGrad);%,'riseTime',obj.dG);
            
            %Append RF and Gradients
            seq.addBlock(GR_seg1,GS_seg1);
            
            seq.addBlock(RF90); % -----------RF 1
            
            seq.addBlock(GR_seg2,GS_seg2);
            
            seq.addBlock(RF90); % -----------RF 2
            
            seq.addBlock(GR_seg3,GS_seg3);
            
            seq.addBlock(RF90); % -----------RF 3
            
            seq.addBlock(GP_seg4);
            
            
        end
        
        
        
    end
    
end