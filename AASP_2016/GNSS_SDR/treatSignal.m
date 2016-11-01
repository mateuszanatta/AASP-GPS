%% treatSignal.m -- starts treatment of GPS data
% Autors:   Alexandre Serio Buscher ==========
%           Hugo Oliveira da Silva ===========
fprintf(['\n\n             SDR GNSS\n',...
        'Initializing Satellites Treatment\n',...
        '----------------------------------------\n'])


%% Generation of noise for each satellite =================================

%fprintf('Noise added\n')
if(~skip)
    fprintf('\t-Noise  . . . . . . . . . . . . . . added;\n');
end
%% Generation of multipath components =====================================
if(~skip)
    fprintf('\t-Multipath components . . . . . . . added;\n');
end
%% Gereneration of receiver signal ========================================
% call genRcvSignal to add all signals and save in approp. file
if(~skip)
    genRcvSignal(satSignal, satellites, settings);
    fprintf('\t-Receiver signal  . . . . . . . . . created;\n');
    %memory cleaning
    clear d; clear ii; clear samplesPerCode; clear satCAtable;
    clear satNAVtable; clear satPtable; %clear satSignal;
    
    disp(['-> Signal Generator is over (elapsed time ', ...
                                        datestr(now - totalTime, 13), ')'])
    clear startTime; clear totalTime;
    fprintf('\n\n');
else
    fprintf('-> Receiver signal  . . . . . . . . . loaded;\n');
end %if(skip)



