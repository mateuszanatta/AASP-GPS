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
freq = set.satL1freq;
%freq = (set.satL1freq + currSat.DoppErr);
%freq = (set.IF + currSat.DoppErr);
carrSin = zeros(1,points);
carrCos = zeros(1,points);
time = 0:1/(set.nyquistGapgen*freq):...
    (0.001*set.nrMSgen-1/(set.nyquistGapgen*freq));

samptime = 0:1:(points-1);
samptime = samptime/(set.nyquistGapgen*set.samplingFreq); %sampling points of generation
sample = 1;

%% Doppler shift rate
%by now no rate, just difference

%% Signal Formation, as GNSS book Appendix B ==============================
fprintf('. |')
ref = 154*set.nyquistGapgen; %Pcode basis = 154 clocks of carrier
TIME = length(time);
marker200us = fix(TIME*0.025);
marker1ms = fix(TIME*0.1);
for now = 1:TIME
    % BPSK itself trough look-up tables
    if( time(now) >= samptime(sample) ) % comparision of generation time with sampling time
        % cossine signal formation
        carrCos(sample) = CAtable(contCA)*NAVtable(contNAV);
        % sine signal formation
        carrSin(sample) = Ptable(contP)*NAVtable(contNAV);
                             % increase sample marker; in the last passage...
        if (sample < points) 
            sample = sample+1; %... it extrapolates index, so I ignored...
        else
            fprintf('|')
            break;
        end %if(sample < points)
        
        %if (~mod(sample,marker200us)) %Ellapsed time marker = "200us"
        %    fprintf('.')
        %end %if(~mod(sample,points/(5*set.nrMSgen)))
     
    end %if( ~mod(now/tt(now)) )
    
    %look-up tables update; all is done by the falling-edge
    if(~mod(now,ref)) % 154 is the number of cicles of P code
       
        if (contP < 10230)
            contP = contP + 1;
            
            if(clockNAV50 < 50)
                clockNAV50 = clockNAV50 + 1;
                
                if(contNAV < 1500)
                    contNAV = contNAV + 1;
                else
                    contNAV = 1;
                end % if(contNAV < 1500)
            
            else
                clockNAV50 = 1;
            end %if(clockNAV50 < 50)
        else
            contP = 1;
   
        end %(contP == 10230)
        
        if(~mod(now,ref*10)) %CA count runs 10 times slower then Pcode
            if (contCA < 1023)
                contCA = contCA + 1;
            else
                contCA = 1;
                fprintf('§')      % markes each complete CA with "§"
            end %(contCA < 1023)     
        end %if(~mod(now,ref*10)
        
    end %if(~mod(now,ref))
    
    
    if(~mod(now,marker200us))
        fprintf('.')
    end
    if(~mod(now,marker1ms))
        fprintf('|')
    end
end % for now = 1:length(time)

%carriers addition

signal = carrCos.*cos(2*pi*freq*samptime);
signal = signal + ...
    carrSin.*sin(2*pi*freq*samptime)*db2mag(-3); %normal attenuation block
fprintf(' .')

end % end function
