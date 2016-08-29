function [mSignal_Delayed] = correlator_banks(pTau, pDelay, pCACode, pPulse, pCodeLength, ...
                            pCodePeriod, pNumberChips, pAmountSatellite, pBankName, pIsConc)
                        
%   Create the correlator bank signals to be used to build the bank of
%   delay correlators. Then, saves it on disk for future use.

%% Input:
%   pTau             - A vector with sampling points for the whole code word of 1023 bits
%   pDelay           - A vector with set of different delays for the Pseudo Random Sequences
%   pCACode          - The Coarse Acquisition Code - Satellite Gold Code words
%   pPulse           - The pulse shape
%   pCodeLength      - Code Length (GPS signals uses a 1023 code length) 
%   pCodePeriod      - The duration of one code period
%   pNumberChips     - Number of chips (duration of pulse)
%   pAmountSatellite - The number of satellites transmitting
%   pBankName        - Define the correlator bank name without '.mat'
%   pIsConc          - Defines whether the signal is concatenated or not.
%                    (Default 0)

%% Output:
%   mSignal_Delayed      - This is the Correlator Bank. An array with a set 
%                        of signals with several pulse delays.

%%   
    % Verify whether the pIsConc paramter was passed
    
    if nargin < 10
        pIsConc = 0;
    end 
    
    % Add a '.mat' at the end of the correlator bank name
    pBankNameFile = strcat(pBankName, '.mat');

%% Generate correlator banks
    
% Firstly we make sure there is no correlator bank saved in the search path
    if exist(pBankNameFile, 'file')
        %If there is a '.mat' file then load it
        
        %When loading the .mat file to a variable, the load method will put
        %the loaded variable into a struct
        mSignal_Delayed = load(pBankNameFile);
        
        %Thus, we need to use getfield to get back the value inside the
        %struct.
        mSignal_Delayed = getfield(mSignal_Delayed, pBankName);
    else
        %if the file cannot be found, then create a new correlator bank
        
        % Initialize the mSignal_Delayed with zeros
        if(pIsConc == 0)
           mSignal_Delayed = zeros(2046, length(pDelay), pAmountSatellite);
        else
           mSignal_Delayed = zeros(22506, length(pDelay), pAmountSatellite);
        end
        
        %Create correlator bank
        for ii = 1:pAmountSatellite;
            for kk = 1:length(pDelay);
                
                if(pIsConc == 0)
                    % Generate correlator bank to signal
                    mSignal_Delayed(:,kk,ii) = gen_signal(pTau + pDelay(kk), ...
                                                          pCACode(:,ii), pPulse, ...
                                                          pCodeLength, pCodePeriod, pNumberChips);
                else
                    % Generate correlator bank to concatenated signals
                    mSignal_Delayed(:,kk,ii) = gen_signal2(pTau + pDelay(kk), ...
                                                          pCACode(:,ii), pPulse, ...
                                                          pCodeLength, pCodePeriod, pNumberChips);
                end
            end
        end
        
        %Defines the variable name to be stored in the .mat file.
        sVarName.(pBankName) = mSignal_Delayed;
        
        %Defines the path to save the correlator bank
        pBankNamePath = strcat('CorrelatorBanks/', pBankNameFile);
        save(pBankNamePath, '-struct', 'sVarName');
    end

end

