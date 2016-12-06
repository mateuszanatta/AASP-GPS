
%% INIT.M
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------
%
%Script initializes settings and environment of the software receiver.
%Then the processing is started.

%--------------------------------------------------------------------------
% CVS record:
% $Id: init.m,v 1.14.2.21 2006/08/22 13:46:00 dpl Exp $

%% Clean up the environment first =========================================
clear; close all; clc;

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions

%% Print startup ==========================================================
fprintf(['\n',...
    'Welcome to:  softGNSS\n\n', ...
    'An open source GNSS SDR software project initiated by:\n\n', ...
    '              Danish GPS Center/Aalborg University\n\n', ...
    'The code was improved by GNSS Laboratory/University of Colorado.\n\n',...
    'The software receiver softGNSS comes with ABSOLUTELY NO WARRANTY;\n',...
    'for details please read license details in the file license.txt. This\n',...
    'is free software, and  you  are  welcome  to  redistribute  it under\n',...
    'the terms described in the license.\n \n',...
    'Adapted by Alexandre Serio Buscher - Instituto Militar de Engenharia\n',...
    '                                   - Fraunhofer EMS - TU Ilmenau\n', ...
    'Sep 2016 \n\n']);
fprintf('                   -------------------------------\n\n');

%% Initialize constants, settings =========================================
settings = initSettings();
fprintf('Struct "settings" initialized\n')
%% Ask for settings structure change ======================================
fprintf('Do you wish to run "setSettings.m"?\n')
choose = input('   1 = Yes    2 = No \n');
fprintf('\n')
if(choose == 1)
    setSettings;
end

%% Generate plot of raw data and ask if ready to start processing =========
choose = 0;
while (choose ~= 1 && choose ~=2 && choose ~= 3)
fprintf('\nChoose data:\n\t')
prompt = (['1 = Dfalt Real Data;  2 = Saved Receiver Data;  ',...
    '3 = Generate Data;\n']);
choose = input(prompt);
    switch(choose)
        case 1      
        case 2 % use generated data
            [settings.fileName,settings.path] = uigetfile('*.bin',...
                'Select the Receiver Data .bin file',...
                ['D:\ALEXANDRE\Vida_Militar_novo\',...
                '2016_IME\Intercambio_Exterior_novo\',...
                'Projetos_Instituto_Fraunhofer_novo\',...
                'Proj_GNSS_AntArray_JoaoPaulo\MatLab_AASP_2016\',...
                'AASP_2016\SatelliteSignals']);
            addpath('satelliteFunctions',settings.path)  % Signal generation functions
            load savedSAT.mat           % SAT matrix of this saved signal
            %load savedSatSignal.mat
            skip = 1; %Skip generation
            generateSignal
        case 3 %generate Data
            settings.path = uigetdir(['D:\ALEXANDRE\Vida_Militar_novo\',...
                '2016_IME\Intercambio_Exterior_novo\',...
                'Projetos_Instituto_Fraunhofer_novo\',...
                'Proj_GNSS_AntArray_JoaoPaulo\MatLab_AASP_2016\',...
                'AASP_2016\SatelliteSignals'],...
                'Select Directory to Save Generated Signal');
            addpath('satelliteFunctions',settings.path)  % Signal generation functions
            settings.fileName = input(...
                'Enter the Receiver Data file saving name:\n','s');
            load SAT.mat                % if you already have a ready SAT matrx
            skip = 0; %Don't skip generation
            generateSignal
        otherwise 
            fprintf('\nBad choice\n');
    end %end while choose
end
fprintf('\n')
try
    fprintf('Probing data (%s)...\n',settings.fileName)
    probeData(settings,choose);
catch
    % There was an error, print it and exit
    errStruct = lasterror;
    disp(errStruct.message);
    %disp('  (run setSettings or change settings in "initSettings.m" to reconfigure)')    
    return;
end %end try    
disp('  Raw IF data loaded and plotted ')
%disp('  (run setSettings or change settings in "initSettings.m" to reconfigure)')
disp(' ');

%%=========================================================================
% change the value of settings.skipTracking to use the pre-correlation
% block or the pre-correlation and post-correlation block

%Utilize this line to process the whole information block of 37ms (original)
%postProcessing

%%=========================================================================
%Utilize this line to process the partial information block of 11ms 
postShortProc

%%=========================================================================