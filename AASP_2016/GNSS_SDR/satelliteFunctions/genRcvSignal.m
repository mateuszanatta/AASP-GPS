function [] = genRcvSignal( signal, sats , settings)
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

%% Antenna Gain ===========================================================
%rcvSignal = db2mag(30)*rcvSignal;

%% Signal Processing ======================================================
% generate front end desired Bandpass Filter (pg 55) but digital
% '' strings; n = cheby2 order    N  low high  fs
n = 8; %for bandpass filter the order is doubled
ftype = 'bandpass';
%Type II Chebyshev Filter
Rs = 30; %attenuation of stopband, between 0 and ripple
Ws = [(settings.IF - settings.BW/2) (settings.IF + settings.BW/2)]...
    /(settings.samplingFreq / 2);
[z,p,k] = cheby2(n,Rs,Ws,ftype);

sos = zp2sos(z,p,k); %numeric precision, [b,a] not good

rcvSignal = sosfilt(sos,rcvSignal);
   
%% Save .bin file for postProcessing script adequated use    
fileNameStr = cat(2,settings.path,settings.fileName,'.bin');
fid = fopen(fileNameStr,'w');
fwrite(fid,rcvSignal,settings.dataType);
fclose(fid);
%save(cat(2,settings.path,'\rcvSignal.mat'),'rcvSignal');

end %function [] = genRcvSignal( signal, sats , set)

