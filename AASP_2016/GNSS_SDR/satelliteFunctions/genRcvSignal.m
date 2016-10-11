function [] = genRcvSignal( signal, sats , set, skip)
%GENRCVSIGNAL 
% INPUTS:   signal        - the satellites generated signal with respectly noise;
%           sats          - satellite struct with their informations
%           set           - settings struct with save file path;

%   add of DoA of LOS and DoA of NLOS, sum of the signals and save as a
%   binary file

%% Addition of DoA ========================================================

%% Received Signal ========================================================
[lin,col] = size(signal);
    rcvSignal = zeros(1,col);
    for ii=1:lin
        rcvSignal = rcvSignal + signal(ii,:);
    end % for ii=1:d
    
if(~skip)    
    
    fileNameStr = cat(2,set.path,set.fileName);
    fid = fopen(fileNameStr,'w');
    fwrite(fid,rcvSignal);
    fclose(fid);
    save(cat(2,set.path,'rcvSignal.mat'),'rcvSignal');
else
    %save('rcvSignal.mat','rcvSignal');
end %function [] = genRcvSignal( signal, sats , set)

