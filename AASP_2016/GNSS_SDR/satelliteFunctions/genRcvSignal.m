function [] = genRcvSignal( signal, sats , set)
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
rcvSignal = db2mag(30)*rcvSignal;

%% Signal Processing ======================================================
% generate front end desired Bandpass Filter (pg 55) but digital
% '' strings; n = cheby2 order    N  low high  fs
n = 20; %for bandpass filter the order is doubled
ftype = 'bandpass';
%{ 
%Type II Chebyshev Filter
Rs = 100; %attenuation of stopband, between 0 and ripple
Ws = [(set.satL1freq - set.BW/2) (set.satL1freq + set.BW/2)]...
    /(set.samplingFreq*42);
[z,p,k] = cheby2(n,Rs,Ws,ftype);
%}

%Type I Chebyshev Filter
Rp = 2;
Wp = [(set.satL1freq - set.BW/2) (set.satL1freq + set.BW/2)]...
    /(set.samplingFreq*set.nyquistGapgen);
[z,p,k] = cheby1(n,Rp,Wp,ftype);

sos = zp2sos(z,p,k);

rcvSignal = sosfilt(sos,rcvSignal);
   
%% Save .bin file for postProcessing script adequated use    
fileNameStr = cat(2,set.path,'\',set.fileName,'.bin');
fid = fopen(fileNameStr,'w');
fwrite(fid,rcvSignal,set.dataType);
fclose(fid);
save(cat(2,set.path,'\rcvSignal.mat'),'rcvSignal');

end %function [] = genRcvSignal( signal, sats , set)

