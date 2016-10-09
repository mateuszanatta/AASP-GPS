function showSatStatus( sats )
%Prints the status of all used Satellites in a table.
%
%showChannelStatus(channel)
%
%   Inputs:
%       channel     - data for each channel. It is used to initialize and
%                   at the processing of the signal (tracking part).
%       settings    - receiver settings

fprintf('\n   *===*=======*===========*===========*=========*======*=======*\n');
fprintf(  ' / |PRN|  DoA  |Code Phase | Dopp. Err.|   SNR   |M.Path|  SIR  |\n');
fprintf(  '*==*===*=======*===========*===========*=========*======*=======*\n');

for satNr = 1 : length(sats)
    fprintf('|%2d| %2d| %3.2fº| %4d chips|  %5d Hz | %2.2fdB|  %2d  | %2.2fdB|\n', ...
                satNr, ...
                sats(satNr).PRN, ...
                sats(satNr).DoA, ...
                sats(satNr).CodPhase, ...
                sats(satNr).DoppErr,...
                sats(satNr).SNR,...
                sats(satNr).Mpath,...
                sats(satNr).SIR);
end

fprintf('*==*===*=======*===========*===========*=========*======*=======*\n\n');

end

