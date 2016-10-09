function [signal] = genSatSignal(currSat, set,...
        points, CAtable, NAVtable, Ptable)
%receiver
% INPUTS:   sat         - is the current satellite;
%           options     - the settings data structure
%           CAtable     - CA look-up table for the satellite
%           NAVtable    - Navigation data look-up table for the satellite
%           Ptable      - P(y) code look-up table for the satellite
% OUTPUTS:  signal      - formed section of signal outputed from this satellite 

%% clocks/counters initial settings =======================================
% the carrier freq is about 1.5 GHz, the fastest clock; but there's no need
% to generate it and since it demands lots of memory; 
fprintf('. ')
%masterCont = 1; % global counter, "the Carrier clock" 
contP = currSat.CodPhase*10; %10.23MHz P Code Clock Counter: 1 inc = 154 masterCont
contCA = currSat.CodPhase; %1.023MHz C/A Code Clock Counter: 1 inc = 10 ContP10_23Mhz
clockNAV50 = 1; %50Hz NavCode Clock: 1 inc = 10230 ContP10_23Mhz
contNAV = 1; % NavMessage Counter                             

%% Carriers generation ====================================================
% Frequency off-set
freq = (set.satL1freq + currSat.DoppErr);
tt = 0:1:(points-1);
tt = tt/(set.nyquistGapgen*set.samplingFreq);
sample = 1;
%CosL1Carr = zeros(1,points);
%SinL1Carr = zeros(1,points);
%CosL1Carr = cos(2*pi*freq*tt/(set.nyquistGapgen*set.samplingFreq));
%SinL1Carr = sin(2*pi*freq*tt/(set.nyquistGapgen*set.samplingFreq));

% Doppler shift rate
%by now no rate, just difference

%% Signal Formation, as GNSS book Appendix B ==============================
fprintf('. ')
signal = zeros(1,points);
time = 0:1/(freq*set.nyquistGapgen):(0.001*set.nrMSgen);
for now = 1:length(time)
    % signal itself trough look-up tables
    if( time(now) >= tt(sample) ) % add um && (sample <= points)
        % cossine signal formation
        signal(sample) = CAtable(contCA)*NAVtable(contNAV)*cos(2*pi*freq*tt(sample));
        % sine signal formation
        signal(sample) = signal(sample) + ...
            ( Ptable(contP)*NAVtable(contNAV)*sin(2*pi*freq*tt(sample)) )*...
            db2mag(-3); %normal attenuation block
        if (sample < points)
            sample = sample+1;
        end %if(sample < points)
        if (~mod(sample,points/10)) % Ellapsed time marker
            fprintf('.')
        end %(~mod(sample,points/20)) % Ellapsed time marker
        
    end %if( ~mod(now/tt(now)) )
    
    %look-up tables update
    if(~mod(now,154)) 
        
        if (contP == 10230)
            contP = 1;
            
            if(clockNAV50 == 50)
                clockNAV50 = 1;
                
                if(contNAV == 1500)
                    contNAV = 1;
                else contNAV = contNAV+1;
                end % if(contNAV > 1499)
            
            else clockNAV50 = clockNAV50+1;
            end %if(clockNAV50 == 50)
        else contP = contP+1;
        end %(contP == 10230)
        
        if(~mod(now,1540))
            if (contCA == 1023)
                contCA = 1;
            else contCA = contCA+1;
            end %(contCA == 1023)
        end %if(~mod(now,1540))
        
        
    end %if(~mod(now,154))
 
end % for now = 1:length(time)
fprintf(' . ')
%% Attenuation - desired SNR ==============================================
signal = signal*db2mag(currSat.SNR);

end % end function
