%% clean up
clc;
close all;
clear;

warning('off','all')
warning

%% code parameters

% code length (GPS: 1023) - GNSS Applications use code lenghts of 1023 bits
Lc=1023;

% duration of one code period (1ms)
Ti=1.0e-3;

% number of chips (chip = bit, in CDMA schemes bits are 1 or -1, called chips) that superimpose (pulse is assumed to be zero after Nc chips)
Nc=4;

% Amount of transmitting simultaneously;
satAmount = 3;

% gps satellite id (1-36) - Identifies satellite signature, after code 36 code words start to repeat itself
% sat_id = round((36 - 1)*rand(1,satAmount) + 1);
sat_id = [7 13 29];

% bandwidth receiver front-end (one-sided) (Transmission Frequency in Hertz 1.026 Mega Hertz
B=1.023e6;

lambda = 3e8/B;
delta = lambda/2;

% oversampling factor pulse (sampling precision of pulse relative to bandwidth)
rho_pulse=100.0;

%% C/D parameters 

% sample time interval
Ts=1/(2*B);

% observed code periods (GPS: 1-20 periods) - For this simulation only one code period is observed - concatenation might apply afterwards
Cp=1;

% observation time
To=Cp*Ti;

%% transmission parameters

% Signal to Noise Ratio (dB)
SNR = -30:2:-20;

% Set of different delays for the Pseudo Random Sequences used for the bank of delay correlators estimation afterwards
shifted_delays_PR = (-0.5:0.25:0.5)*Ts;

% Set of different delays used for differente Line Of Sight (LOS) possible Signals - These delays are used to compare signals and filter the transmitted ones 
LOS_delays = (-1:0.01:1)*Ts;

% Pulse delay used for the Transmitted Signal - Most relevant transmission parameter
% delaySet = (1+1)*rand(1,satAmount) - 1;
% delaySet = round(delaySet*100)/100;
% signal_tx_delay = delaySet*Ts;
signal_tx_delay = [-0.73 0.27 0.81]*Ts;


%% generate signals

% generate C/A code (Coarse Acquisition Code - Satellite Gold Code words)
id = zeros(1023, satAmount);
for ii = 1:satAmount;
    id(:,ii) = gen_ca_code(sat_id(ii));
end

% generate pulse shape (rectangular pulse after band-pass filter) - Nyquits pulse tha modulates one chip (bit)
pulse = gen_pulse_rect_bp(B,rho_pulse,Lc,Ti,Nc);

% sample points - Sampling points for the whole code word of 1023 bits
tau = (0.0:Ts:To-Ts)';

% Set of Signals with different pulse delays within the same Pseudo Random sequence used for a satellite transmission
signal_shiftedPR = zeros(2046, length(shifted_delays_PR), satAmount);
for ii = 1:satAmount;
    for kk = 1:length(shifted_delays_PR);
        % generate correlators bank signals of Pseudo Random Sequences
        signal_shiftedPR(:,kk,ii) = gen_signal(tau + shifted_delays_PR(kk), id(:,ii), pulse, Lc, Ti, Nc);
    end
end

signal_shiftedPR2 = zeros(22506, length(shifted_delays_PR), satAmount);
for ii = 1:satAmount;
    for kk = 1:length(shifted_delays_PR);
        % generate correlators bank signals of Pseudo Random Sequences
        signal_shiftedPR2(:,kk,ii) = gen_signal2(tau + shifted_delays_PR(kk), id(:,ii), pulse, Lc, Ti, Nc);
    end
end

% Set of possible Line of Sight signals, each one with a different pulse
% delay - for further estimation
signal_LOS_copies = zeros(2046, length(shifted_delays_PR), satAmount);
for ii = 1:satAmount;
    for kk = 1:length(LOS_delays)
        signal_LOS_copies(:,kk,ii) = gen_signal(tau + LOS_delays(kk), id(:,ii), pulse, Lc, Ti, Nc);
    end
end

signal_LOS_copies2 = zeros(22506, length(shifted_delays_PR), satAmount);
for ii = 1:satAmount;
    for kk = 1:length(LOS_delays)
        signal_LOS_copies2(:,kk,ii) = gen_signal2(tau + LOS_delays(kk), id(:,ii), pulse, Lc, Ti, Nc);
    end
end

phaseTx = [-63 -2 71]*pi/180;
phaseTx = sort(phaseTx, 'descend');
signalPhases = exp(1j*phaseTx);

% Base for DOA RMSE calculus - variable to help calculation afterwards
rmse_temp_phase = repmat(phaseTx, length(SNR), 1);

% Transmitted Signal
argumentTx = zeros(2046, satAmount);
for ii = 1:satAmount;
    argumentTx(:,ii) = gen_signal(tau+signal_tx_delay(ii),id(:,ii),pulse,Lc,Ti,Nc)*signalPhases(ii);
end

argumentTx2 = zeros(22506, satAmount);
for ii = 1:satAmount;
    argumentTx2(:,ii) = gen_signal2(tau+signal_tx_delay(ii),id(:,ii),pulse,Lc,Ti,Nc)*signalPhases(ii);
end

% For an unique antenna model signals will be transmitted through an additive channel
signal_tx = sum(argumentTx,2);

% Concatenation For an unique antenna model signals will be transmitted through an additive channel
signal_tx_concat = sum(argumentTx2,2);

% Multiple Antennas Signal Model (d x S, matrix)
signal_tx_MA = argumentTx;

% Multiple Antennas Signal Model (d x S, matrix)
signal_tx_MA_concat = argumentTx2;

%% Multiple Antennas Scheme
% Amount of Antennas in the array
M1 = 4;
M2 = 16;

% Generate Antena Array model - Steering Matrix
A1 = array_matr_Rd(phaseTx, M1, 1);
A2 = array_matr_Rd(phaseTx, M2, 1);

% Receveid Signal with multiple antennas - Noise still to be added
X0_1 = A1*transpose(signal_tx_MA);
X0_2 = A2*transpose(signal_tx_MA);
X0_conc = A2*transpose(signal_tx_MA_concat);
X0_conc_5 = A1*transpose(signal_tx_MA_concat);


%% RMSE
iterationsRMSE = 1; % Amount of Iterations for RMSE calculation (Root Mean Square Error - determines algorithm performance throughout executions)
for idrmse = 1:iterationsRMSE;
    % Important to no put ";" in here to check how fast the simulation is running
    idrmse
    for idxSNR = 1:length(SNR);
        %% AWGN transmission

        % White Noise standard deviation - SNR calculation
        sigma = 10^(-SNR(idxSNR)/20)/sqrt(2);

        % Draw ZMCSCG (zero mean) WHITE noise with variance sigma^2
        whiteNoise = (randn(size(signal_tx))+ 1*1i*randn(size(signal_tx)))*sigma;

        % Received Signal
        signal_rx = signal_tx + whiteNoise;
        
        % Draw ZMCSCG (zero mean) WHITE noise with variance sigma^2 - CONC
        whiteNoise_conc = (randn(size(signal_tx_concat))+ 1*1i*randn(size(signal_tx_concat)))*sigma;
        % Received Signal one antenna CONC
        signal_rx_conc = signal_tx_concat + whiteNoise_conc;
        
        % Matricial Noise
        whiteNoise_M1 = (randn(size(X0_1))+ 1*1i*randn(size(X0_1)))*sigma;
        
        whiteNoise_M2 = (randn(size(X0_2))+ 1*1i*randn(size(X0_2)))*sigma;
        
        whiteNoise_conc = (randn(size(X0_conc))+ 1*1i*randn(size(X0_conc)))*sigma;
        
        whiteNoise_conc_5 = (randn(size(X0_conc_5))+ 1*1i*randn(size(X0_conc_5)))*sigma;
        
        % Received Signal Multiple Antennas
        % 5
        X_5 = X0_1 + whiteNoise_M1;
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 CONC
        X_conc_5 = X0_conc_5 + whiteNoise_conc_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16
        X_16 = X0_2 + whiteNoise_M2;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 CONC
        X_conc_16 = X0_conc + whiteNoise_conc;
                
        % Estimating the noisy subspace - Signal Processing through SVD
        % 5
        Us_5 = est_sigsubsp_classic(X_5,satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 CONC
        Us_conc_5 = est_sigsubsp_classic(X_conc_5,satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16
        Us_2 = est_sigsubsp_classic(X_16,satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 CONC
        Us_conc_16 = est_sigsubsp_classic(X_conc_16,satAmount);
       
        % DOA (Direction of Arrival) Estimation
        % 5
        phase_est_5 = standard_esprit_Rd(Us_5,M1);
        phase_est_5 = sort(phase_est_5,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 conc
        phase_est_conc_5 = standard_esprit_Rd(Us_conc_5,M1);
        phase_est_conc_5 = sort(phase_est_conc_5,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 16
        phase_est_16 = standard_esprit_Rd(Us_2,M2);
        phase_est_16 = sort(phase_est_16,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 Conc
        phase_est_conc_16 = standard_esprit_Rd(Us_conc_16,M2);
        phase_est_conc_16 = sort(phase_est_conc_16,'descend');
        
        % Estimated DOAS for multiple antennas with different SNRs
        % ----------------------------------------------------------------------------------------------------------------------
        phase_MA_est_5(idxSNR, :) = phase_est_5;
        phase_MA_est_conc_5(idxSNR, :) = phase_est_conc_5;
        
        % ----------------------------------------------------------------------------------------------------------------------
        phase_MA_est_16(idxSNR, :) = phase_est_16;
        phase_MA_est_conc_16(idxSNR, :) = phase_est_conc_16;
        
        % Once the phi_est is computed we can generate back the matriz A
       % ----------------------------------------------------------------------------------------------------------------------
        A_est_5 = array_matr_Rd(phase_est_5,M1);
        A_est_conc_5 = array_matr_Rd(phase_est_conc_5,M1);
        
        for ii = 1:satAmount;
            T = GetFocusedVIT(phase_est_5(ii)*180/pi, delta, M1, 0.8);
            X_vit_5 = T*X_5;
            U_vit_5 = est_sigsubsp_classic(X_vit_5, 1);
            phase_vit_5(ii) = standard_esprit_Rd(U_vit_5,M1);
            a_vit_5 = array_matr_Rd(phase_vit_5(ii),M1);
            S_vit_5(ii,:) = a_vit_5'*X_vit_5;
            
%             [ Pw, theta, PwMusic ] = calculaCapon( X_5, M1, delta, lambda );
%             figure;
%             plot(theta, Pw, 'r', theta, PwMusic);
%             
%             [ Pw2, theta2, PwMusic2 ] = calculaCapon( X_vit_5, M1, delta, lambda );
%             figure;
%             plot(theta2, Pw2,'r', theta2, PwMusic2);
        end
        
        % ----------------------------------------------------------------------------------------------------------------------
        A_est_16 = array_matr_Rd(phase_est_16,M2);
        A_est_conc_16 = array_matr_Rd(phase_est_conc_16,M2);
        
        for ii = 1:satAmount;
            T = GetFocusedVIT(phase_est_16(ii)*180/pi, 0.5, M2, 0.8);
            X_vit_16 = T*X_16;
            
            a_vit_16 = array_matr_Rd(0,M2);
            S_vit_16(ii,:) = a_vit_16'*X_vit_16;
        end
        
        % Generate Signal Matrix - Moore-Penrose Pseudo Inverse
        % ----------------------------------------------------------------------------------------------------------------------
        S_est_5 = pinv(A_est_5)*X_5;
        
        %Generate 4 antenna Signal Matrix with the already known Steering
        %Matrix
        S_5 = pinv(A1)*X_5;
         
        S_est_conc_5 = pinv(A_est_conc_5)*X_conc_5;
        
        % ----------------------------------------------------------------------------------------------------------------------
        S_est_16 = pinv(A_est_16)*X_16;
        
        %Generate 16 antenna Signal Matrix with the already known Steering
        %Matrix
        S_16 = pinv(A2)*X_16;
        
        S_est_conc_16 = pinv(A_est_conc_16)*X_conc_16;

        %% SPS/FBA Design
        
        % Apply Forward Backward Averaging to Signal with multiple Antennas
        
        % 5 SPS + FBA + Concat
        X_fba_5 = dofba_tensor(X_conc_5);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA
        X_fba_5 = dofba_tensor(X_5);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS + FBA + CONC
        X_fba = dofba_tensor(X_conc_16);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        X_fba_16 = dofba_tensor(X_16);
        

        % L determines how the antennas in the array are selected
        L1 = 1;
        L2 = 1;
        
        % Apply Spatial Smoothing through Signal
        % 5 SPS + FBA + CONC
        X_fba_SS_5 = spatialsmooth_meastensor(X_fba_5,L1);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS
        X_SS_5 = spatialsmooth_meastensor(X_5,L1);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS + FBA + CONC
        X_fba_SS = spatialsmooth_meastensor(X_fba,L2);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS
        X_SS_16 = spatialsmooth_meastensor(X_16,L2);
        
        % Low Rank approximation through SVD
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA SPS CONC
        Us_fba_SS_5 = est_sigsubsp_classic(X_fba_SS_5, satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA 
        Us_fba_5 = est_sigsubsp_classic(X_fba_5, satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS 
        Us_SS_5 = est_sigsubsp_classic(X_SS_5, satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA SPS CONC
        Us_fba_SS = est_sigsubsp_classic(X_fba_SS, satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        Us_fba_16 = est_sigsubsp_classic(X_fba_16, satAmount);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS
        Us_SS_16 = est_sigsubsp_classic(X_SS_16, satAmount);

        % Estimated Direction of Arrivals
        % 5 SPS + FBA CONC
        phi_fba_SS_est_5 = standard_esprit_Rd(Us_fba_SS_5,M1-L1+1);
        phi_fba_SS_est_5 = sort(phi_fba_SS_est_5,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA
        phi_fba_5 = standard_esprit_Rd(Us_fba_5,M1-L1+1);
        phi_fba_5 = sort(phi_fba_5,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS
        phi_SS_5 = standard_esprit_Rd(Us_SS_5,M1-L1+1);
        phi_SS_5 = sort(phi_SS_5,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA SPS CONC
        phi_fba_SS_est_16 = standard_esprit_Rd(Us_fba_SS,M2-L2+1);
        phi_fba_SS_est_16 = sort(phi_fba_SS_est_16,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        phi_fba_16 = standard_esprit_Rd(Us_fba_16,M2-L2+1);
        phi_fba_16 = sort(phi_fba_16,'descend');
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS
        phi_SS_16 = standard_esprit_Rd(Us_SS_16,M2-L2+1);
        phi_SS_16 = sort(phi_SS_16,'descend');
        
        
        % Estimated DOAS for multiple antennas with different SNRs (SPS/FBA)
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS + FBA + CONC
        phase_MA_fba_SS_est_5(idxSNR,:) = phi_fba_SS_est_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA
        phase_MA_fba_5(idxSNR,:) = phi_fba_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS
        phase_MA_SS_5(idxSNR,:) = phi_SS_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA SPS CONC
        phase_MA_fba_SS_est(idxSNR,:) = phi_fba_SS_est_16;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        phase_MA_fba_16(idxSNR,:) = phi_fba_16;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS
        phase_MA_SS_16(idxSNR,:) = phi_SS_16;
        
        %         for ii = 1:satAmount;
            
            % colocar a doa
            T = GetFocusedVIT(phi_fba_SS_est_16(2), 0.5, M2, 0.8);

            X_vit_16 = T*X_conc_16;

%             Us_vit = est_sigsubsp_classic(X_vit,satAmount);

            % DOA (Direction of Arrival) Estimation
%             phase_est_vit = standard_esprit_Rd(Us_vit,M);
%             phase_est_vit(ii) = min(abs(phase_est_vit));
%             phase_est_vit_corrected(ii) = phase_est_vit(ii) + phase_est(ii);
%         end
        
%         phase_MA_fba_SS_est(idxSNR,:) = phase_est_vit_corrected;

        % Once the phi_est is computed we can generate back the matriz A
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA SPS CONC
        A_est_fba_SS_5 = array_matr_Rd(phi_fba_SS_est_5,M1);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA
        A_est_fba_5 = array_matr_Rd(phi_fba_5,M1);
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS
        A_est_SS_5 = array_matr_Rd(phi_SS_5,M1);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA SPS CONC
        A_est_fba_SS_16 = array_matr_Rd(phi_fba_SS_est_16,M2);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        A_est_fba_16 = array_matr_Rd(phi_fba_16,M2);
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SPS
        A_est_SS_16 = array_matr_Rd(phi_SS_16,M2);

        % Generate Signal Matrix - Moore-Penrose Pseudo Inverse
        
        % 5 FBA SPS CONC
        S_est_fba_SS_5 = pinv(A_est_fba_SS_5)*X_conc_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 FBA
        S_est_fba_5 = pinv(A_est_fba_5)*X_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 SPS
        S_est_SS_5 = pinv(A_est_SS_5)*X_5;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA SPS CONC
        S_est_fba_SS_16 = pinv(A_est_fba_SS_16)*X_conc_16;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 FBA
        S_est_fba_16 = pinv(A_est_fba_16)*X_16;
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 SS
        S_est_SS_16 = pinv(A_est_SS_16)*X_16;
        
S_est_fba_SS_vit = pinv(A_est_fba_SS_16)*X_vit_16;
        
%% /////////////             /////////////////
%                DLL Begin??
% /////////////             /////////////////
    %% Lowest Correlation Estimation Noise Free Environment

        % Offset correlators - Unique Antenna Model;
        offAngle = (-90:10:90)*pi/180;
        offSetCorr = exp(-1j*offAngle);

        for ii = 1:satAmount;
            % Correlation Bank Matched to a set of different possible LOS signals
            corr_matched_bank(:,:,ii) = signal_LOS_copies(:,:,ii)'*signal_shiftedPR(:,:,ii);
        end
        
        for ii = 1:satAmount;
            % Correlation Bank Matched to a set of different possible LOS signals
            corr_matched_bank2(:,:,ii) = signal_LOS_copies2(:,:,ii)'*signal_shiftedPR2(:,:,ii);
        end
        [rowArgConc, ~, ~] = size(corr_matched_bank2);

        % Correlation between Transmitted signal and PR signals delayed
        % One tries to find here the lowest value on deviance. This value
        % indicates the most likelihood to the transmitted signal.
        for ii = 1:satAmount;
            for kk = 1:length(offAngle);
                corr_signal_tx_bank(:,:,ii) = (signal_tx*offSetCorr(kk))'*signal_shiftedPR(:,:,ii);

                % Delay Deviance between copies of Bank and Transmitted Signal
                ml_argument = bsxfun(@minus, corr_matched_bank(:,:,ii), corr_signal_tx_bank(:,:,ii));
                deviance_LOS(:,kk,ii) = sum(sqrt(ml_argument.^2),2);  
            end
        end

        % Over this step one will estimate the delay that corresponds to
        % the transmitted one. The lowest value indicates the estimated
        % value.
        for ii = 1:satAmount;
            offSetEstimator = sum(sqrt(imag(deviance_LOS(:,:,ii)).^2));
            [~, offSetIdx] = min(offSetEstimator);

            deviance_LOS_est(:,:,ii) = deviance_LOS(:, offSetIdx,ii);
            phase_estimation_n_free(:,ii) = offAngle(offSetIdx);
        end

        %% Correlator Estimator AWGN Environment
        %%%%%% Same estimation as the one above, but in a Noisy channel
        %%%%%% (AWGN)

        % Correlation between Transmitted signal and PR signals delayed
        [rowArg, ~, ~] = size(corr_matched_bank);
        for ii = 1:satAmount;
            for kk = 1:length(offAngle);
                corr_signal_rx_bank(:,:,ii) = (signal_rx*offSetCorr(kk))'*signal_shiftedPR(:,:,ii);

                % Delay Deviance between copies of Bank and Transmitted Signal
                temporarySub = repmat(corr_signal_rx_bank(:,:,ii), rowArg, 1);
                ml_argument_rx = corr_matched_bank(:,:,ii) - temporarySub;
                deviance_LOS_rx(:,kk,ii) = sum(sqrt(ml_argument_rx.^2),2);  
            end
        end
        
        % Over this step one will estimate the delay that corresponds to
        % the transmitted one. The lowest value indicates the estimated
        % value.
        for ii = 1:satAmount;
            offSetEstimator = sum(sqrt(imag(deviance_LOS_rx(:,:,ii)).^2));
            [~, offSetIdx] = min(offSetEstimator);

            deviance_LOS_est_rx(:,:,ii,idxSNR) = deviance_LOS_rx(:, offSetIdx,ii);
            
            % Estimating transmitted signal phases. This algorithm does not
            % consist in a classic DOA estimation algorithm, but only a Maximum
            % Likelihood Estimator to find the phase as well.
            phase_estimation(:,ii) = offAngle(offSetIdx);
            phase_est_DOA(:,ii,idxSNR) = offAngle(offSetIdx);
        end
        
        phaseTx_est(idxSNR,:) = phase_estimation;
        
 %---------------------------------------------------------------------------------
 % One antenna - CONC
        
        % Correlation between Transmitted signal and PR signals delayed
        [rowArg, ~, ~] = size(corr_matched_bank);
        for ii = 1:satAmount;
            for kk = 1:length(offAngle);
                corr_signal_rx_bank_conc(:,:,ii) = (signal_rx_conc*offSetCorr(kk))'*signal_shiftedPR2(:,:,ii);

                % Delay Deviance between copies of Bank and Transmitted Signal
                temporarySub_conc = repmat(corr_signal_rx_bank_conc(:,:,ii), rowArg, 1);
                ml_argument_rx_conc = corr_matched_bank(:,:,ii) - temporarySub_conc;
                deviance_LOS_rx_conc(:,kk,ii) = sum(sqrt(ml_argument_rx_conc.^2),2);  
            end
        end
        
        % Over this step one will estimate the delay that corresponds to
        % the transmitted one. The lowest value indicates the estimated
        % value.
        for ii = 1:satAmount;
            offSetEstimator_conc = sum(sqrt(imag(deviance_LOS_rx_conc(:,:,ii)).^2));
            [~, offSetIdx_conc] = min(offSetEstimator_conc);

            deviance_LOS_est_rx_conc(:,:,ii,idxSNR) = deviance_LOS_rx_conc(:, offSetIdx_conc,ii);
            
            % Estimating transmitted signal phases. This algorithm does not
            % consist in a classic DOA estimation algorithm, but only a Maximum
            % Likelihood Estimator to find the phase as well.
            phase_estimation_conc(:,ii) = offAngle(offSetIdx_conc);
            phase_est_DOA_conc(:,ii,idxSNR) = offAngle(offSetIdx_conc);
        end
        
        phaseTx_est_conc(idxSNR,:) = phase_estimation_conc;
       
        %% Correlator Esimator Multiple Antennas       
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS
        % Correlation between Estimated signal and PR signals delayed   
        S_calc_1 = transpose(S_est_5);
        for ii = 1:satAmount;
        	corr_S_est_bank_1(:,:,ii) = (S_calc_1(:,ii)*exp(-1j*phase_est_5(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_MA_1 = repmat(corr_S_est_bank_1(:,:,ii), rowArg, 1);
            S_est_argument_rx_1 = corr_matched_bank(:,:,ii) - tempSub_MA_1;
            deviance_LOS_S_est_1(:,:,ii) = sum(sqrt(S_est_argument_rx_1.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_1 = sum(sqrt(imag(deviance_LOS_S_est_1(:,:,ii)).^2));
            [~, offSetIdx_S_est_1] = min(offSetEstimator_S_est_1);

            deviance_LOS_S_est_MA_5(:,:,ii,idxSNR) = deviance_LOS_S_est_1(:, offSetIdx_S_est_1,ii);
        end
        % 5 ANTENAS - Known Steering Matrix
        % Correlation between Estimated signal and PR signals delayed Using the already Knwon Steering Matrix   
        S_Known = transpose(S_5);
        for ii = 1:satAmount;
        	corr_S_Known_bank_1(:,:,ii) = (S_Known(:,ii)*exp(-1j*phase_est_5(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_Known_MA_5 = repmat(corr_S_Known_bank_1(:,:,ii), rowArg, 1);
            S_Known_argument_rx_5 = corr_matched_bank(:,:,ii) - tempSub_Known_MA_5;
            deviance_LOS_S_Known_5(:,:,ii) = sum(sqrt(S_Known_argument_rx_5.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_Known_5 = sum(sqrt(imag(deviance_LOS_S_Known_5(:,:,ii)).^2));
            [~, offSetIdx_S_Known_5] = min(offSetEstimator_S_Known_5);

            deviance_LOS_S_Known_MA_5(:,:,ii,idxSNR) = deviance_LOS_S_Known_5(:, offSetIdx_S_Known_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS CONC
        % Correlation between Estimated signal and PR signals delayed   
        S_calc_conc_5 = transpose(S_est_conc_5);
        for ii = 1:satAmount;
        	corr_S_est_bank_conc_5(:,:,ii) = (S_calc_conc_5(:,ii)*exp(-1j*phase_est_conc_5(ii)))'*signal_shiftedPR2(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_MA_conc_5 = repmat(corr_S_est_bank_conc_5(:,:,ii), rowArgConc, 1);
            S_est_argument_rx_conc_5 = corr_matched_bank2(:,:,ii) - tempSub_MA_conc_5;
            deviance_LOS_S_est_conc_5(:,:,ii) = sum(sqrt(S_est_argument_rx_conc_5.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_conc_5 = sum(sqrt(imag(deviance_LOS_S_est_conc_5(:,:,ii)).^2));
            [~, offSetIdx_S_est_conc_5] = min(offSetEstimator_S_est_conc_5);

            deviance_LOS_S_est_MA_conc_5(:,:,ii,idxSNR) = deviance_LOS_S_est_conc_5(:, offSetIdx_S_est_conc_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS VIT
        % Correlation between Estimated signal and PR signals delayed   
        S_calc_vit_5 = transpose(S_vit_5);
        for ii = 1:satAmount;
        	corr_S_est_bank_vit_5(:,:,ii) = (S_calc_vit_5(:,ii)*exp(-1j*phase_vit_5(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_MA_vit_5 = repmat(corr_S_est_bank_vit_5(:,:,ii), rowArg, 1);
            S_est_argument_rx_vit_5 = corr_matched_bank(:,:,ii) - tempSub_MA_vit_5;
            deviance_LOS_S_est_vit_5(:,:,ii) = sum(sqrt(S_est_argument_rx_vit_5.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_vit_5 = sum(sqrt(imag(deviance_LOS_S_est_vit_5(:,:,ii)).^2));
            [~, offSetIdx_S_est_vit_5] = min(offSetEstimator_S_est_vit_5);

            deviance_LOS_S_est_MA_vit_5(:,:,ii,idxSNR) = deviance_LOS_S_est_vit_5(:, offSetIdx_S_est_vit_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 ANTENAS
        S_calc_16 = transpose(S_est_16);
        for ii = 1:satAmount;
        	corr_S_est_bank_16(:,:,ii) = (S_calc_16(:,ii)*exp(-1j*phase_est_16(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_MA_16 = repmat(corr_S_est_bank_16(:,:,ii), rowArg, 1);
            S_est_argument_rx_16 = corr_matched_bank(:,:,ii) - tempSub_MA_16;
            deviance_LOS_S_est_16(:,:,ii) = sum(sqrt(S_est_argument_rx_16.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_16 = sum(sqrt(imag(deviance_LOS_S_est_16(:,:,ii)).^2));
            [~, offSetIdx_S_est_16] = min(offSetEstimator_S_est_16);

            deviance_LOS_S_est_MA_16(:,:,ii,idxSNR) = deviance_LOS_S_est_16(:, offSetIdx_S_est_16,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 ANTENAS - Known Steering Matrix
        S_Known_16 = transpose(S_16);
        for ii = 1:satAmount;
        	corr_S_Known_bank_16(:,:,ii) = (S_Known_16(:,ii)*exp(-1j*phase_est_16(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_Known_MA_16 = repmat(corr_S_Known_bank_16(:,:,ii), rowArg, 1);
            S_Known_argument_rx_16 = corr_matched_bank(:,:,ii) - tempSub_Known_MA_16;
            deviance_LOS_S_Known_16(:,:,ii) = sum(sqrt(S_Known_argument_rx_16.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_Known_16 = sum(sqrt(imag(deviance_LOS_S_Known_16(:,:,ii)).^2));
            [~, offSetIdx_S_Known_16] = min(offSetEstimator_S_Known_16);

            deviance_LOS_S_Known_MA_16(:,:,ii,idxSNR) = deviance_LOS_S_Known_16(:, offSetIdx_S_Known_16,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 ANTENAS CONC
        S_calc_conc = transpose(S_est_conc_16);
        for ii = 1:satAmount;
        	corr_S_est_bank_conc(:,:,ii) = (S_calc_conc(:,ii)*exp(-1j*phase_est_conc_16(ii)))'*signal_shiftedPR2(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            tempSub_MA_conc = repmat(corr_S_est_bank_conc(:,:,ii), rowArgConc, 1);
            S_est_argument_rx_conc = corr_matched_bank2(:,:,ii) - tempSub_MA_conc;
            deviance_LOS_S_est_conc(:,:,ii) = sum(sqrt(S_est_argument_rx_conc.^2),2);
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_conc = sum(sqrt(imag(deviance_LOS_S_est_conc(:,:,ii)).^2));
            [~, offSetIdx_S_est_conc] = min(offSetEstimator_S_est_conc);

            deviance_LOS_S_est_MA_conc(:,:,ii,idxSNR) = deviance_LOS_S_est_conc(:, offSetIdx_S_est_conc,ii);
        end
        
        %% Correlator Esimator Multiple Antennas FBA + SPS
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS FBA SPS CONC
        % Correlation between Estimated signal and PR signals delayed
        S_calc_fba_SS_5 = transpose(S_est_fba_SS_5);
        for ii = 1:satAmount;
            corr_S_est_fba_SS_bank_5(:,:,ii) = (S_calc_fba_SS_5(:,ii)*exp(-1j*phi_fba_SS_est_5(ii)))'*signal_shiftedPR2(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_fba_SS_argument_rx_5 = bsxfun(@minus, corr_matched_bank2(:,:,ii), corr_S_est_fba_SS_bank_5(:,:,ii));
            deviance_LOS_S_est_fba_SS_5(:,:,ii) = sum(sqrt(S_est_fba_SS_argument_rx_5.^2),2);  
            
        end
        
        % Over this step one will estimate the delay that corresponds to
        % the transmitted one. The lowest value indicates the estimated
        % value.
        for ii = 1:satAmount;
            offSetEstimator_S_est_fba_SS_5 = sum(sqrt(imag(deviance_LOS_S_est_fba_SS_5(:,:,ii)).^2));
            [~, offSetIdx_S_est_fba_SS_5] = min(offSetEstimator_S_est_fba_SS_5);

            deviance_LOS_S_est_fba_SS_MA_5(:,:,ii,idxSNR) = deviance_LOS_S_est_fba_SS_5(:, offSetIdx_S_est_fba_SS_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS FBA
        % Correlation between Estimated signal and PR signals delayed
        S_calc_fba_5 = transpose(S_est_fba_5);
        for ii = 1:satAmount;
            corr_S_est_fba_bank_5(:,:,ii) = (S_calc_fba_5(:,ii)*exp(-1j*phi_fba_5(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_fba_argument_rx_5 = bsxfun(@minus, corr_matched_bank(:,:,ii), corr_S_est_fba_bank_5(:,:,ii));
            deviance_LOS_S_est_fba_5(:,:,ii) = sum(sqrt(S_est_fba_argument_rx_5.^2),2);  
            
        end
       
        for ii = 1:satAmount;
            offSetEstimator_S_est_fba_5 = sum(sqrt(imag(deviance_LOS_S_est_fba_5(:,:,ii)).^2));
            [~, offSetIdx_S_est_fba_5] = min(offSetEstimator_S_est_fba_5);

            deviance_LOS_S_est_fba_MA_5(:,:,ii,idxSNR) = deviance_LOS_S_est_fba_5(:, offSetIdx_S_est_fba_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 5 ANTENAS SS
        % Correlation between Estimated signal and PR signals delayed
        S_calc_SS_5 = transpose(S_est_SS_5);
        for ii = 1:satAmount;
            corr_S_est_SS_bank_5(:,:,ii) = (S_calc_SS_5(:,ii)*exp(-1j*phi_SS_5(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_SS_argument_rx_5 = bsxfun(@minus, corr_matched_bank(:,:,ii), corr_S_est_SS_bank_5(:,:,ii));
            deviance_LOS_S_est_SS_5(:,:,ii) = sum(sqrt(S_est_SS_argument_rx_5.^2),2);  
            
        end
       
        for ii = 1:satAmount;
            offSetEstimator_S_est_SS_5 = sum(sqrt(imag(deviance_LOS_S_est_SS_5(:,:,ii)).^2));
            [~, offSetIdx_S_est_SS_5] = min(offSetEstimator_S_est_SS_5);

            deviance_LOS_S_est_SS_MA_5(:,:,ii,idxSNR) = deviance_LOS_S_est_SS_5(:, offSetIdx_S_est_SS_5,ii);
        end
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 Antenas FBA SPS CONC
        % Correlation between Estimated signal and PR signals delayed
        S_calc_fba_SS_16 = transpose(S_est_fba_SS_16);
        for ii = 1:satAmount;
            corr_S_est_fba_SS_bank_16(:,:,ii) = (S_calc_fba_SS_16(:,ii)*exp(-1j*phi_fba_SS_est_16(ii)))'*signal_shiftedPR2(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_fba_SS_argument_rx_16 = bsxfun(@minus, corr_matched_bank2(:,:,ii), corr_S_est_fba_SS_bank_16(:,:,ii));
            deviance_LOS_S_est_fba_SS_16(:,:,ii) = sum(sqrt(S_est_fba_SS_argument_rx_16.^2),2);  
            
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_fba_SS_16 = sum(sqrt(imag(deviance_LOS_S_est_fba_SS_16(:,:,ii)).^2));
            [~, offSetIdx_S_est_fba_SS_16] = min(offSetEstimator_S_est_fba_SS_16);

            deviance_LOS_S_est_fba_SS_MA_16(:,:,ii,idxSNR) = deviance_LOS_S_est_fba_SS_16(:, offSetIdx_S_est_fba_SS_16,ii);
        end
        
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 Antenas FBA
        % Correlation between Estimated signal and PR signals delayed
        S_calc_fba_16 = transpose(S_est_fba_16);
        for ii = 1:satAmount;
            corr_S_est_fba_16(:,:,ii) = (S_calc_fba_16(:,ii)*exp(-1j*phi_fba_16(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_fba_argument_rx_16 = bsxfun(@minus, corr_matched_bank(:,:,ii), corr_S_est_fba_16(:,:,ii));
            deviance_LOS_S_est_fba_16(:,:,ii) = sum(sqrt(S_est_fba_argument_rx_16.^2),2);  
            
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_fba_16 = sum(sqrt(imag(deviance_LOS_S_est_fba_16(:,:,ii)).^2));
            [~, offSetIdx_S_est_fba_16] = min(offSetEstimator_S_est_fba_16);

            deviance_LOS_S_est_fba_MA_16(:,:,ii,idxSNR) = deviance_LOS_S_est_fba_16(:, offSetIdx_S_est_fba_16,ii);
        end
        
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 Antenas SPS
        % Correlation between Estimated signal and PR signals delayed
        S_calc_SS_16 = transpose(S_est_SS_16);
        for ii = 1:satAmount;
            corr_S_est_SS_16(:,:,ii) = (S_calc_SS_16(:,ii)*exp(-1j*phi_SS_16(ii)))'*signal_shiftedPR(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_SS_argument_rx_16 = bsxfun(@minus, corr_matched_bank(:,:,ii), corr_S_est_SS_16(:,:,ii));
            deviance_LOS_S_est_SS_16(:,:,ii) = sum(sqrt(S_est_SS_argument_rx_16.^2),2);  
            
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_SS_16 = sum(sqrt(imag(deviance_LOS_S_est_SS_16(:,:,ii)).^2));
            [~, offSetIdx_S_est_SS_16] = min(offSetEstimator_S_est_SS_16);

            deviance_LOS_S_est_SS_MA_16(:,:,ii,idxSNR) = deviance_LOS_S_est_SS_16(:, offSetIdx_S_est_fba_16,ii);
        end
        
        % ----------------------------------------------------------------------------------------------------------------------
        % 16 Antenas VIT
        % Correlation between Estimated signal and PR signals delayed
        S_calc_fba_SS_vit = transpose(S_est_fba_SS_vit);
        for ii = 1:satAmount;
            corr_S_est_fba_SS_bank_vit(:,:,ii) = (S_calc_fba_SS_vit(:,ii)*exp(-1j*phi_fba_SS_est_16(ii)))'*signal_shiftedPR2(:,:,ii);

            % Delay Deviance between copies of Bank and Transmitted Signal
            S_est_fba_SS_argument_rx_vit = bsxfun(@minus, corr_matched_bank2(:,:,ii), corr_S_est_fba_SS_bank_vit(:,:,ii));
            deviance_LOS_S_est_fba_SS_vit(:,:,ii) = sum(sqrt(S_est_fba_SS_argument_rx_vit.^2),2);  
            
        end
        
        for ii = 1:satAmount;
            offSetEstimator_S_est_fba_SS_vit = sum(sqrt(imag(deviance_LOS_S_est_fba_SS_vit(:,:,ii)).^2));
            [~, offSetIdx_S_est_fba_SS_vit] = min(offSetEstimator_S_est_fba_SS_vit);

            deviance_LOS_S_est_fba_SS_MA_vit(:,:,ii,idxSNR) = deviance_LOS_S_est_fba_SS_vit(:, offSetIdx_S_est_fba_SS_vit,ii);
        end
    end
    %% /////////////             /////////////////
    %                DLL End??
    % /////////////             /////////////////
    %%  DOA RMSE
    
    % Root Mean Square Error for Unique Antenna Model
    rmse_DOA_est_1(:,:,idrmse) =  (phaseTx_est - rmse_temp_phase).^2;
    
    % Root Mean Square Error for Unique Antenna Model - CONC
    rmse_DOA_est_1_conc(:,:,idrmse) =  (phaseTx_est_conc - rmse_temp_phase).^2;
    
    % Root Mean Square Error for Multiple Antenna Model
    rmse_DOA_MA_est_5(:,:,idrmse) =  (phase_MA_est_5 - rmse_temp_phase).^2;
    rmse_DOA_MA_est_fba_5(:,:,idrmse) =  (phase_MA_fba_5 - rmse_temp_phase).^2;
    rmse_DOA_MA_est_SS_5(:,:,idrmse) =  (phase_MA_SS_5 - rmse_temp_phase).^2;
    rmse_DOA_MA_est_conc_5(:,:,idrmse) =  (phase_MA_est_conc_5 - rmse_temp_phase).^2;
    rmse_DOA_MA_fba_SS_est_5(:,:,idrmse) =  (phase_MA_fba_SS_est_5 - rmse_temp_phase).^2;
    
    rmse_DOA_MA_est_16(:,:,idrmse) =  (phase_MA_est_16 - rmse_temp_phase).^2;
    rmse_DOA_MA_fba_16(:,:,idrmse) =  (phase_MA_fba_16 - rmse_temp_phase).^2;
    rmse_DOA_MA_SS_16(:,:,idrmse) =  (phase_MA_SS_16 - rmse_temp_phase).^2;
    rmse_DOA_MA_est_conc_16(:,:,idrmse) =  (phase_MA_est_conc_16 - rmse_temp_phase).^2;
    rmse_DOA_MA_fba_SS_est_16(:,:,idrmse) =  (phase_MA_fba_SS_est - rmse_temp_phase).^2;
    
    % Root Mean Square Error for Multiple Antenna Model plus SPS and FBA
    rmse_DOA_MA_fba_SS_est_vit(:,:,idrmse) =  (phase_MA_fba_SS_est - rmse_temp_phase).^2;
        
    %% Delay Estimation No Scheme per SNR

    % Find indexes of pulse delay to estimate signals with different SNRs
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_est] = min(real(deviance_LOS_est_rx(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_est(ii,snr) = LOS_delays(indx_est); 
        end
    end
    
%     delay_est = roundn(delay_est/Ts,-1)*Ts;
    
    % Calculate RMSE inner argument
    argMin = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg(:,:,ii,snr,idrmse) = (real(delay_est) - argMin).^2;
        end
    end
    
% Find indexes of pulse delay to estimate signals with different SNRs -
% CONC
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_est_conc] = min(real(deviance_LOS_est_rx_conc(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_est_conc(ii,snr) = LOS_delays(indx_est_conc); 
        end
    end
    
        % Calculate RMSE inner argument - CONC
    argMin_conc = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_conc(:,:,ii,snr,idrmse) = (real(delay_est_conc) - argMin_conc).^2;
        end
    end
    
    %% Delay Estimation Multiple Antennas
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 ANTENAS
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_5] = min(real(deviance_LOS_S_est_MA_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_est_5(ii,snr) = LOS_delays(indx_S_est_5); 
        end
    end

    
    % Calculate RMSE inner argument
    argMin_S_est_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_5(:,:,ii,snr,idrmse) = (real(delay_S_est_5) - argMin_S_est_5).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 ANTENAS - Known Steering Matrix
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_Known_5] = min(real(deviance_LOS_S_Known_MA_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_Known_5(ii,snr) = LOS_delays(indx_S_Known_5); 
        end
    end

    
    % Calculate RMSE inner argument
    argMin_S_Known_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_Known_5(:,:,ii,snr,idrmse) = (real(delay_S_Known_5) - argMin_S_Known_5).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 ANTENAS VIT
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_vit_5] = min(real(deviance_LOS_S_est_MA_vit_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_est_vit_5(ii,snr) = LOS_delays(indx_S_est_vit_5); 
        end
    end

    
    % Calculate RMSE inner argument
    argMin_S_est_vit_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_vit_5(:,:,ii,snr,idrmse) = (real(delay_S_est_vit_5) - argMin_S_est_vit_5).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 ANTENAS CONC
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_conc_5] = min(real(deviance_LOS_S_est_MA_conc_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_est_conc_5(ii,snr) = LOS_delays(indx_S_est_conc_5); 
        end
    end

    
    % Calculate RMSE inner argument
    argMin_S_est_conc_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_conc_5(:,:,ii,snr,idrmse) = (real(delay_S_est_conc_5) - argMin_S_est_conc_5).^2;
        end
    end
    % ----------------------------------------------------------------------------------------------------------------------
    % 16 ANTENAS
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_16] = min(real(deviance_LOS_S_est_MA_16(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_est_16(ii,snr) = LOS_delays(indx_S_est_16); 
        end
    end
    
    
    % Calculate RMSE inner argument
    argMin_S_est_16 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_16(:,:,ii,snr,idrmse) = (real(delay_S_est_16) - argMin_S_est_16).^2;
        end
    end
    % ----------------------------------------------------------------------------------------------------------------------
    % 16 ANTENAS - Known Steering Matrix
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_Known_16] = min(real(deviance_LOS_S_Known_MA_16(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_Known_16(ii,snr) = LOS_delays(indx_S_Known_16); 
        end
    end
    
    
    % Calculate RMSE inner argument
    argMin_S_Known_16 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_Known_16(:,:,ii,snr,idrmse) = (real(delay_S_Known_16) - argMin_S_Known_16).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    % 16 ANTENAS CONC
    % Find delays for Multiple Antenna Schemes
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_conc] = min(real(deviance_LOS_S_est_MA_conc(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_est_conc(ii,snr) = LOS_delays(indx_S_est_conc); 
        end
    end

    
    % Calculate RMSE inner argument
    argMin_S_est_conc = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_conc(:,:,ii,snr,idrmse) = (real(delay_S_est_conc) - argMin_S_est_conc).^2;
        end
    end
    
    %% Delay Estimation Multiple Antennas FBA + SS
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 ANTENAS FBA + SS + CONC
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_1_5] = min(real(deviance_LOS_S_est_fba_SS_MA_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_fba_SS_est_5(ii,snr) = LOS_delays(indx_S_est_1_5); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_fba_SS_est_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_fba_SS_5(:,:,ii,snr,idrmse) = (real(delay_S_fba_SS_est_5) - argMin_S_fba_SS_est_5).^2;
        end
    end
    % ----------------------------------------------------------------------------------------------------------------------
    % 5 Antenna Schemes FBA
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_fba_5] = min(real(deviance_LOS_S_est_fba_MA_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_fba_est_5(ii,snr) = LOS_delays(indx_S_est_fba_5); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_fba_est_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_fba_5(:,:,ii,snr,idrmse) = (real(delay_S_fba_est_5) - argMin_S_fba_est_5).^2;
        end
    end
    
     % ----------------------------------------------------------------------------------------------------------------------
    % 5 Antenna Schemes SPS
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_SS_5] = min(real(deviance_LOS_S_est_SS_MA_5(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_SS_est_5(ii,snr) = LOS_delays(indx_S_est_SS_5); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_SS_est_5 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_SS_5(:,:,ii,snr,idrmse) = (real(delay_S_SS_est_5) - argMin_S_SS_est_5).^2;
        end
    end
    % ----------------------------------------------------------------------------------------------------------------------
    %  16 FBA + SS + CONC
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_SPS_FBA_16] = min(real(deviance_LOS_S_est_fba_SS_MA_16(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_fba_SS_est_16(ii,snr) = LOS_delays(indx_S_est_SPS_FBA_16); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_fba_SS_est_16 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_fba_SS_16(:,:,ii,snr,idrmse) = (real(delay_S_fba_SS_est_16) - argMin_S_fba_SS_est_16).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    %  16 FBA
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_fba_16] = min(real(deviance_LOS_S_est_fba_MA_16(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_fba_est_16(ii,snr) = LOS_delays(indx_S_est_fba_16); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_fba_est_16 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_fba_16(:,:,ii,snr,idrmse) = (real(delay_S_fba_est_16) - argMin_S_fba_est_16).^2;
        end
    end
    % ----------------------------------------------------------------------------------------------------------------------
    %  16 SPS
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_SS_16] = min(real(deviance_LOS_S_est_SS_MA_16(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_SS_est_16(ii,snr) = LOS_delays(indx_S_est_SS_16); 
        end
    end
    
    % Calculate RMSE inner argument
    argMin_S_SS_est_16 = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_SS_16(:,:,ii,snr,idrmse) = (real(delay_S_SS_est_16) - argMin_S_SS_est_16).^2;
        end
    end
    
    % ----------------------------------------------------------------------------------------------------------------------
    % 16 antenas VIT
    % Find delays for Multiple Antenna Schemes FBA + SS
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            [~, indx_S_est_SPS_FBA_vit] = min(real(deviance_LOS_S_est_fba_SS_MA_vit(:,:,ii,snr)));
            
            % Delays retrieved here allow transmitted signal
            % reconstruction for different SNRs
            delay_S_fba_SS_est_vit(ii,snr) = LOS_delays(indx_S_est_SPS_FBA_vit); 
        end
    end

%     delay_S__fba_SS_est = roundn(delay_S_est/Ts, -1)*Ts;
    
    % Calculate RMSE inner argument
    argMin_S_fba_SS_est_vit = repmat(signal_tx_delay', 1, length(SNR));
    for snr = 1:length(SNR);
        for ii = 1:satAmount;
            rmseArg_S_est_fba_SS_vit(:,:,ii,snr,idrmse) = (real(delay_S_fba_SS_est_vit) - argMin_S_fba_SS_est_vit).^2;
        end
    end
end
% Loop END

%% Delays RMSE Calculation

% Unique Antenna Delay RMSE calculation
RMSE = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp = rmseArg(:,:,ii,snr,idrmse);
            RMSE = RMSE + argTemp;
        end
    end
end
RMSE = sqrt(RMSE/iterationsRMSE);
% -------------------------------------------------------------------------------------------------------
% Unique Antenna Delay RMSE calculation - Conc
RMSE_conc = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_conc = rmseArg_conc(:,:,ii,snr,idrmse);
            RMSE_conc = RMSE_conc + argTemp_conc;
        end
    end
end
RMSE_conc = sqrt(RMSE_conc/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS
RMSE_S_est_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_5 = rmseArg_S_est_5(:,:,ii,snr,idrmse);
            RMSE_S_est_5 = RMSE_S_est_5 + argTemp_S_est_5;
        end
    end
end
RMSE_S_est_5 = sqrt(RMSE_S_est_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS - Known Steering
% Matrix
RMSE_S_Known_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_Known_5 = rmseArg_S_Known_5(:,:,ii,snr,idrmse);
            RMSE_S_Known_5 = RMSE_S_Known_5 + argTemp_S_Known_5;
        end
    end
end
RMSE_S_Known_5 = sqrt(RMSE_S_Known_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS VIT
RMSE_S_est_vit_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_vit_5 = rmseArg_S_est_vit_5(:,:,ii,snr,idrmse);
            RMSE_S_est_vit_5 = RMSE_S_est_vit_5 + argTemp_S_est_vit_5;
        end
    end
end
RMSE_S_est_vit_5 = sqrt(RMSE_S_est_vit_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS CONC
RMSE_S_est_conc_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_conc_5 = rmseArg_S_est_conc_5(:,:,ii,snr,idrmse);
            RMSE_S_est_conc_5 = RMSE_S_est_conc_5 + argTemp_S_est_conc_5;
        end
    end
end
RMSE_S_est_conc_5 = sqrt(RMSE_S_est_conc_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS FBA
RMSE_S_est_fba_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_fba_5 = rmseArg_S_est_fba_5(:,:,ii,snr,idrmse);
            RMSE_S_est_fba_5 = RMSE_S_est_fba_5 + argTemp_S_est_fba_5;
        end
    end
end
RMSE_S_est_fba_5 = sqrt(RMSE_S_est_fba_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 5 ANTENAS SPS
RMSE_S_est_SS_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_SS_5 = rmseArg_S_est_SS_5(:,:,ii,snr,idrmse);
            RMSE_S_est_SS_5 = RMSE_S_est_SS_5 + argTemp_S_est_SS_5;
        end
    end
end
RMSE_S_est_SS_5 = sqrt(RMSE_S_est_SS_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation with FBA and SPS 5 ANTENAS
RMSE_S_est_fba_SS_5 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_fba_SS_5 = rmseArg_S_est_fba_SS_5(:,:,ii,snr,idrmse);
            RMSE_S_est_fba_SS_5 = RMSE_S_est_fba_SS_5 + argTemp_S_est_fba_SS_5;
        end
    end
end
RMSE_S_est_fba_SS_5 = sqrt(RMSE_S_est_fba_SS_5/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 16 ANTENAS
RMSE_S_est_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_16 = rmseArg_S_est_16(:,:,ii,snr,idrmse);
            RMSE_S_est_16 = RMSE_S_est_16 + argTemp_S_est_16;
        end
    end
end
RMSE_S_est_16 = sqrt(RMSE_S_est_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 16 ANTENAS - Known Steering
% Matrix
RMSE_S_Known_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_Known_16 = rmseArg_S_Known_16(:,:,ii,snr,idrmse);
            RMSE_S_Known_16 = RMSE_S_Known_16 + argTemp_S_Known_16;
        end
    end
end
RMSE_S_Known_16 = sqrt(RMSE_S_Known_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 16 ANTENAS CONC
RMSE_S_est_conc_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_conc = rmseArg_S_est_conc(:,:,ii,snr,idrmse);
            RMSE_S_est_conc_16 = RMSE_S_est_conc_16 + argTemp_S_est_conc;
        end
    end
end
RMSE_S_est_conc_16 = sqrt(RMSE_S_est_conc_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 16 ANTENAS FBA
RMSE_S_est_fba_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_fba_16 = rmseArg_S_est_fba_16(:,:,ii,snr,idrmse);
            RMSE_S_est_fba_16 = RMSE_S_est_fba_16 + argTemp_S_est_fba_16;
        end
    end
end
RMSE_S_est_fba_16 = sqrt(RMSE_S_est_fba_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation 16 ANTENAS SPS
RMSE_S_est_SS_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_SS_16 = rmseArg_S_est_SS_16(:,:,ii,snr,idrmse);
            RMSE_S_est_SS_16 = RMSE_S_est_SS_16 + argTemp_S_est_SS_16;
        end
    end
end
RMSE_S_est_SS_16 = sqrt(RMSE_S_est_SS_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation with FBA and SPS 16 antenas
RMSE_S_est_fba_SS_16 = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_fba_SS = rmseArg_S_est_fba_SS_16(:,:,ii,snr,idrmse);
            RMSE_S_est_fba_SS_16 = RMSE_S_est_fba_SS_16 + argTemp_S_est_fba_SS;
        end
    end
end
RMSE_S_est_fba_SS_16 = sqrt(RMSE_S_est_fba_SS_16/iterationsRMSE);
% ----------------------------------------------------------------------------------------------------------------------
% Multiple Antennas Delay RMSE calculation with FBA and SPS 16 antenas vit
RMSE_S_est_fba_SS_vit = zeros(satAmount,length(SNR));
for ii = 1:satAmount;
    for snr = 1:length(SNR);
        for idrmse = 1:iterationsRMSE;
            argTemp_S_est_fba_SS_vit = rmseArg_S_est_fba_SS_vit(:,:,ii,snr,idrmse);
            RMSE_S_est_fba_SS_vit = RMSE_S_est_fba_SS_vit + argTemp_S_est_fba_SS_vit;
        end
    end
end
RMSE_S_est_fba_SS_vit = sqrt(RMSE_S_est_fba_SS_vit/iterationsRMSE);

%% DOA RMSE
RMSE_DOA = sqrt(sum(rmse_DOA_est_1,3)/iterationsRMSE);
RMSE_DOA_conc = sqrt(sum(rmse_DOA_est_1_conc,3)/iterationsRMSE);
RMSE_DOA_MA_5 = sqrt(sum(rmse_DOA_MA_est_5,3)/iterationsRMSE);
RMSE_DOA_MA_conc_5 = sqrt(sum(rmse_DOA_MA_est_conc_5,3)/iterationsRMSE);
RMSE_DOA_MA_fba_5 = sqrt(sum(rmse_DOA_MA_est_fba_5,3)/iterationsRMSE);
RMSE_DOA_MA_SS_5 = sqrt(sum(rmse_DOA_MA_est_SS_5,3)/iterationsRMSE);
RMSE_DOA_MA_FBA_SS_5 = sqrt(sum(rmse_DOA_MA_fba_SS_est_5,3)/iterationsRMSE);


RMSE_DOA_MA_16 = sqrt(sum(rmse_DOA_MA_est_16,3)/iterationsRMSE);
RMSE_DOA_MA_conc_16 = sqrt(sum(rmse_DOA_MA_est_conc_16,3)/iterationsRMSE);
RMSE_DOA_MA_fba_16 = sqrt(sum(rmse_DOA_MA_fba_16,3)/iterationsRMSE);
RMSE_DOA_MA_SS_16 = sqrt(sum(rmse_DOA_MA_SS_16,3)/iterationsRMSE);
RMSE_DOA_MA_FBA_SS_16 = sqrt(sum(rmse_DOA_MA_fba_SS_est_16,3)/iterationsRMSE);
RMSE_DOA_MA_FBA_SS_VIT = sqrt(sum(rmse_DOA_MA_fba_SS_est_vit,3)/iterationsRMSE);

%% DELAY ALL PLOTS
%Set plot style
set(0,'DefaultTextFontName','Times New Roman',...
'DefaultTextFontSize',16,...
'DefaultAxesFontName','Times New Roman',...
'DefaultAxesFontSize',16,...
'DefaultLineLineWidth',3.5,...
'DefaultLineMarkerSize',10.75);
figure();
mean_RMSE = sum(RMSE)/3;
semilogy(SNR, mean_RMSE,'o-', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 0]);
hold on;

mean_RMSE_conc = sum(RMSE_conc)/3;
semilogy(SNR, mean_RMSE_conc,'+--', 'LineWidth', 2, ...
'MarkerFaceColor', [0.5 0.5 0.5]);

mean_RMSE_S_est_5 = sum(RMSE_S_est_5)/3;
semilogy(SNR, mean_RMSE_S_est_5,'*-.', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 0.5]);

mean_RMSE_Known_est_5 = sum(RMSE_S_Known_5)/3;
semilogy(SNR, mean_RMSE_Known_est_5,'o-.', 'LineWidth', 2, ...
'MarkerFaceColor', [1 1 0.5]);

mean_RMSE_S_est_SS_5 = sum(RMSE_S_est_SS_5)/3;
semilogy(SNR, mean_RMSE_S_est_SS_5,'x:', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0.5 0]);

mean_RMSE_S_est_fba_5 = sum(RMSE_S_est_fba_5)/3;
semilogy(SNR, mean_RMSE_S_est_fba_5,'s-', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0.5 0.5]);
%     semilogy(SNR, sum(RMSE_S_est_vit_5)/3,'o-.', 'LineWidth', 2, ...
%                 'MarkerFaceColor', [0.5 0 0]);

mean_RMSE_S_est_conc_5 = sum(RMSE_S_est_conc_5)/3;
semilogy(SNR, mean_RMSE_S_est_conc_5,'d--', 'LineWidth', 2, ...
'MarkerFaceColor', [0.5 0 0.5]);

mean_RMSE_S_est_fba_SS_5 = sum(RMSE_S_est_fba_SS_5)/3;
semilogy(SNR, mean_RMSE_S_est_fba_SS_5,'p:', 'LineWidth', 2, ...
'MarkerFaceColor', [0.5 0.5 0]);

mean_RMSE_S_est_16 = sum(RMSE_S_est_16)/3;
semilogy(SNR, mean_RMSE_S_est_16,'h-.', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 1]);

mean_RMSE_S_Known_16 = sum(RMSE_S_Known_16)/3;
semilogy(SNR, mean_RMSE_S_Known_16,'o-.', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 1]);

mean_RMSE_S_est_SS_16 = sum(RMSE_S_est_SS_16)/3;
semilogy(SNR, mean_RMSE_S_est_SS_16,'v--', 'LineWidth', 2, ...
'MarkerFaceColor', [0.2 0.05 0]);

mean_RMSE_S_est_fba_16 = sum(RMSE_S_est_fba_16)/3;
semilogy(SNR, mean_RMSE_S_est_fba_16,'<:', 'LineWidth', 2, ...
'MarkerFaceColor', [0 1 1]);

mean_RMSE_S_est_conc_16 = sum(RMSE_S_est_conc_16)/3;
semilogy(SNR, mean_RMSE_S_est_conc_16,'.--', 'LineWidth', 2, ...
'MarkerFaceColor', [1 0 1]);

mean_RMSE_S_est_fba_SS_16 = sum(RMSE_S_est_fba_SS_16)/3;
semilogy(SNR, mean_RMSE_S_est_fba_SS_16,'>:', 'LineWidth', 2, ...
'MarkerFaceColor', [1 1 0]);
hold off;
grid on;
xlabel('Signal to Noise Ratio (SNR).');
ylabel('RMSE of \tau');
legend('Single antenna',...
'Single antenna - CONC',...
[num2str(M1) ' ULA'],...
[num2str(M1) ' ULA - A Known'],...
[num2str(M1) ' ULA - SPS'],...
[num2str(M1) ' ULA - FBA'],...
[num2str(M1) ' ULA - CONC'],...
[num2str(M1) ' ULA - FBA+SPS+CONC'],...
[num2str(M2) ' ULA'],...
[num2str(M2) ' ULA - A Known'],...
[num2str(M2) ' ULA - SPS'],...
[num2str(M2) ' ULA - FBA'],...
[num2str(M2) ' ULA - CONC'],...
[num2str(M2) ' ULA - FBA+SPS+CONC']);
%         [num2str(M1) ' ULA com VIT'],...
%% DOA ALL PLOTS
figure;

mean_RMSE_DOA = sum(RMSE_DOA')/3;
semilogy(SNR, mean_RMSE_DOA,'o:', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 0]);
hold on;

mean_RMSE_DOA_conc = sum(RMSE_DOA_conc')/3;
semilogy(SNR, mean_RMSE_DOA_conc,'+-', 'LineWidth', 2, ...
'MarkerFaceColor', [0.5 0.5 0.5]);

mean_RMSE_DOA_MA_5 = sum(RMSE_DOA_MA_5')/3;
semilogy(SNR, mean_RMSE_DOA_MA_5,'*--', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0 0.5]);

mean_RMSE_DOA_MA_conc_5 = sum(RMSE_DOA_MA_conc_5')/3;
semilogy(SNR, mean_RMSE_DOA_MA_conc_5,'x-.', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0.5 0]);

mean_RMSE_DOA_MA_fba_5 = sum(RMSE_DOA_MA_fba_5')/3;
semilogy(SNR, mean_RMSE_DOA_MA_fba_5,'s--', 'LineWidth', 2, ...
'MarkerFaceColor', [0 0.5 0.5]);

mean_RMSE_DOA_MA_SS_5 = sum(RMSE_DOA_MA_SS_5')/3;
semilogy(SNR, mean_RMSE_DOA_MA_SS_5,'p-.', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0.5 0 0]);

mean_RMSE_DOA_MA_FBA_SS_5 = sum(RMSE_DOA_MA_FBA_SS_5')/3;
semilogy(SNR, mean_RMSE_DOA_MA_FBA_SS_5,'h:', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0.5 0 0.5]);
% semilogy(SNR, sum(RMSE_DOA')/3,'o:', 'LineWidth', 2, ...
%         'MarkerFaceColor', [0 0 1]);
% hold on;

mean_RMSE_DOA_MA_16 = sum(RMSE_DOA_MA_16')/3;
semilogy(SNR, mean_RMSE_DOA_MA_16,'v-', 'LineWidth', 2, ...
'MarkerFaceColor', [ 1 0 0]);

mean_RMSE_DOA_MA_conc_16 = sum(RMSE_DOA_MA_conc_16')/3;
semilogy(SNR, mean_RMSE_DOA_MA_conc_16,'d:', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0 1 0]);

mean_RMSE_DOA_MA_fba_16 = sum(RMSE_DOA_MA_fba_16')/3;
semilogy(SNR, mean_RMSE_DOA_MA_fba_16,'<--', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0 1 0]);

mean_RMSE_DOA_MA_SS_16 = sum(RMSE_DOA_MA_SS_16')/3;
semilogy(SNR, mean_RMSE_DOA_MA_SS_16,'.:', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0 1 0]);

mean_RMSE_DOA_MA_FBA_SS_16 = sum(RMSE_DOA_MA_FBA_SS_16')/3;
semilogy(SNR, mean_RMSE_DOA_MA_FBA_SS_16,'>:', 'LineWidth', 2, ...
'MarkerFaceColor', [ 0 1 0]);
hold off;
grid on;
xlabel('Signal to Noise Ratio (SNR).');
ylabel('RMSE of \theta');
legend('Single antenna',...
'Single antenna - CONC',...
[num2str(M1) ' ULA'],...
[num2str(M1) ' ULA - SPS'],...
[num2str(M1) ' ULA - FBA'],...
[num2str(M1) ' ULA - CONC'],...
[num2str(M1) ' ULA - FBA+SPS+CONC'],...
[num2str(M2) ' ULA'],...
[num2str(M2) ' ULA - SPS'],...
[num2str(M2) ' ULA - FBA'],...
[num2str(M2) ' ULA - CONC'],...
[num2str(M2) ' ULA - FBA+SPS+CONC']);
% 
% %% RMSE PLOT 5 Antenas
% figure;
% for ii = 1:4;
%     subplot(2,2,ii);
%     if ii ~= 4;
%         semilogy(SNR, RMSE(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 0]);
%         hold on;
%         semilogy(SNR, RMSE_S_est_5(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
% 
%         semilogy(SNR, RMSE_S_est_SS_5(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 0]);
% 
%         semilogy(SNR, RMSE_S_est_fba_5(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 1]);
% 
%     %     semilogy(SNR, RMSE_S_est_vit_5(ii,:),'o-.', 'LineWidth', 2, ...
%     %                 'MarkerFaceColor', [1 0 0]);
% 
%         semilogy(SNR, RMSE_S_est_conc_5(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 0 1]);
% 
%         semilogy(SNR, RMSE_S_est_fba_SS_5(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 1 0]); 
%         hold off;
%     else
%         semilogy(SNR, sum(RMSE)/3,'o-.', 'LineWidth', 2, ...
%                 'MarkerFaceColor', [0 0 0]);
%         hold on;
%         semilogy(SNR, sum(RMSE_S_est_5)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
% 
%         semilogy(SNR, sum(RMSE_S_est_SS_5)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 0]);
% 
%         semilogy(SNR, sum(RMSE_S_est_fba_5)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 1]);
% 
%     %     semilogy(SNR, sum(RMSE_S_est_vit_5)/3,'o-.', 'LineWidth', 2, ...
%     %                 'MarkerFaceColor', [1 0 0]);
% 
%         semilogy(SNR, sum(RMSE_S_est_conc_5)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 0 1]);
% 
%         semilogy(SNR, sum(RMSE_S_est_fba_SS_5)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 1 0]); 
%         hold off;
%     end
%     grid on;
%     xlabel('Signal to Noise Ratio (SNR).');
%     ylabel('RMSE - Delay Error.');
%     legend('No Scheme',...
%         [num2str(M1) ' ULA'],...
%         [num2str(M1) ' ULA - SPS'],...
%         [num2str(M1) ' ULA - FBA'],...
%         [num2str(M1) ' ULA - CONC'],...
%         [num2str(M1) ' ULA - FBA+SPS+CONC']);
% %         [num2str(M1) ' ULA com VIT'],...
%     
% end
% 
% %% RMSE PLOT 16 Antenas
% figure;
% for ii = 1:4;
%     subplot(2,2,ii);
%     if ii ~= 4;
%         semilogy(SNR, RMSE(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 0]);
%         hold on;
%         semilogy(SNR, RMSE_S_est_16(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
% 
%         semilogy(SNR, RMSE_S_est_SS_16(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 0]);
% 
%         semilogy(SNR, RMSE_S_est_fba_16(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 1]);
% 
%         semilogy(SNR, RMSE_S_est_conc_16(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 0 1]);
% 
%         semilogy(SNR, RMSE_S_est_fba_SS_16(ii,:),'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 1 0]); 
%         hold off;
%         grid on;
%     else
%         semilogy(SNR, sum(RMSE)/3,'o-.', 'LineWidth', 2, ...
%                 'MarkerFaceColor', [0 0 0]);
%         hold on;
%         semilogy(SNR, sum(RMSE_S_est_16)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
% 
%         semilogy(SNR, sum(RMSE_S_est_SS_16)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 0]);
% 
%         semilogy(SNR, sum(RMSE_S_est_fba_16)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 1 1]);
% 
%         semilogy(SNR, sum(RMSE_S_est_conc_16)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 0 1]);
% 
%         semilogy(SNR, sum(RMSE_S_est_fba_SS_16)/3,'o-.', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [1 1 0]); 
%         hold off;
%         grid on;
%     end
%     xlabel('Signal to Noise Ratio (SNR).');
%     ylabel('RMSE - Delay Error.');
%     legend('No Scheme', ...
%         [num2str(M2) ' ULA'],...
%         [num2str(M2) ' ULA - SPS'],...
%         [num2str(M2) ' ULA - FBA'],...
%         [num2str(M2) ' ULA - CONC'],...
%         [num2str(M2) ' ULA - FBA+SPS+CONC']);
%     
% end
% 
% %% RMSE DOA PLOT 5 - Antenas
% figure;
% for ii = 1:4;
%     if ii ~= 4;
%         subplot(2,2,ii);
%         semilogy(SNR, RMSE_DOA(:,ii),'o:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
%         hold on;
%         semilogy(SNR, RMSE_DOA_MA_5(:,ii),'s:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 1 0 0]);
%                 
%         semilogy(SNR, RMSE_DOA_MA_conc_5(:,ii),'d:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 0 1 0]);
%                 
%         semilogy(SNR, RMSE_DOA_MA_fba_5(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, RMSE_DOA_MA_SS_5(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, RMSE_DOA_MA_FBA_SS_5(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%     else
%         subplot(2,2,ii);
%         semilogy(SNR, sum(RMSE_DOA')/3,'o:', 'LineWidth', 2, ...
%                 'MarkerFaceColor', [0 0 1]);
%         hold on;
%         semilogy(SNR, sum(RMSE_DOA_MA_5')/3,'s:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 1 0 0]);
%                 
%         semilogy(SNR, sum(RMSE_DOA_MA_conc_5')/3,'d:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 0 1 0]);
%                 
%         semilogy(SNR, sum(RMSE_DOA_MA_fba_5')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, sum(RMSE_DOA_MA_SS_5')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, sum(RMSE_DOA_MA_FBA_SS_5')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%     end
%     grid on;
%     xlabel('Signal to Noise Ratio (SNR).');
%     ylabel('RMSE - DOA Error');
%     legend('No Scheme',...
%         [num2str(M1) ' ULA'],...
%         [num2str(M1) ' ULA - CONC'],...
%         [num2str(M1) ' ULA - FBA'], ...
%         [num2str(M1) ' ULA - SPS'],...
%         [num2str(M1) ' ULA - FBA,SPS e CONC']);
%     hold off;
% end
% 
% %% RMSE DOA PLOT 16 - Antenas
% figure;
% for ii = 1:4;
%     if ii ~= 4;
%         subplot(2,2,ii);
%         semilogy(SNR, RMSE_DOA(:,ii),'o:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [0 0 1]);
%         hold on;
%         semilogy(SNR, RMSE_DOA_MA_16(:,ii),'s:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 1 0 0]);
%                 
%         semilogy(SNR, RMSE_DOA_MA_conc_16(:,ii),'d:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 0 1 0]);
%                 
%         semilogy(SNR, RMSE_DOA_MA_fba_16(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, RMSE_DOA_MA_SS_16(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, RMSE_DOA_MA_FBA_SS_16(:,ii),'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%     else
%         subplot(2,2,ii);
%         semilogy(SNR, sum(RMSE_DOA')/3,'o:', 'LineWidth', 2, ...
%                 'MarkerFaceColor', [0 0 1]);
%         hold on;
%         semilogy(SNR, sum(RMSE_DOA_MA_16')/3,'s:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 1 0 0]);
%                 
%         semilogy(SNR, sum(RMSE_DOA_MA_conc_16')/3,'d:', 'LineWidth', 2, ...
%                     'MarkerFaceColor', [ 0 1 0]);
%                 
%         semilogy(SNR, sum(RMSE_DOA_MA_fba_16')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, sum(RMSE_DOA_MA_SS_16')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%         
%         semilogy(SNR, sum(RMSE_DOA_MA_FBA_SS_16')/3,'d:', 'LineWidth', 2, ...
%             'MarkerFaceColor', [ 0 1 0]);
%     end
%     grid on;
%     xlabel('Signal to Noise Ratio (SNR).');
%     ylabel('RMSE - DOA Error');
%     legend('No Scheme', ...
%         [num2str(M2) ' ULA'],...
%         [num2str(M2) ' ULA - CONC'],...
%         [num2str(M2) ' ULA - FBA'], ...
%         [num2str(M2) ' ULA - SPS'],...
%         [num2str(M2) ' ULA - FBA,SPS e CONC']);
%     hold off;
% end
% 
% 
% 
% %% Delay Plot
% figure;
% for ii = 1:satAmount;
%     plot(LOS_delays/Ts,real(deviance_LOS_est(:,:,ii)), '--', 'LineWidth', 2);
%     hold on;
% end
% hold off;
% grid on;
% title('Maximum Likelihood Estimator');
% legend(['Satellite Number ' num2str(sat_id(1)) ' - Delay Pulse ' num2str(signal_tx_delay(1)/Ts)],...
%        ['Satellite Number ' num2str(sat_id(2)) ' - Delay Pulse ' num2str(signal_tx_delay(2)/Ts)],...
%        ['Satellite Number ' num2str(sat_id(3)) ' - Delay Pulse ' num2str(signal_tx_delay(3)/Ts)]);
% xlabel('Set of Possible Delays \tau (xT_s)');
% ylabel('Deviance of Transmitted Pulse');
% 
% %% Delay Plot AWGN - No Scheme
% figure;
% for ii = [1 2 3];
%     subplot(2,2,ii);
%     plot(LOS_delays/Ts,real(deviance_LOS_est(:,:,ii)));
%     hold on;
%     for snr = [1 6 11 16];
%         plot(LOS_delays/Ts,real(deviance_LOS_est_rx(:,:,ii,snr)));
%     end
%     hold off;
%     grid on;
%     legend('Noise Free Environment', '-30 dB SNR',...
%             '-20 dB SNR','-10 dB SNR', '0 dB SNR');
%     xlabel(['Set of Possible Delays \tau (xT_s) - Signal Delay = ' num2str(signal_tx_delay(ii)/Ts) ' x T_s']);
%     ylabel('Deviance of Transmitted Pulse');
% end
% 
% %% Delay Plot AWGN - 5 Antenas Scheme
% figure;
% for ii = 1:satAmount;
%     subplot(2,2,ii);
%     plot(LOS_delays/Ts,real(deviance_LOS_est(:,:,ii)));
%     hold on;
%     for snr = [1 6 11 16];
%         plot(LOS_delays/Ts,real(deviance_LOS_S_est_MA_5(:,:,ii,snr)));
%     end
%     hold off;
%     grid on;
%     legend('Noise Free Environment', '-30 dB SNR',...
%             '-20 dB SNR', '-10 dB SNR', '0 dB SNR');
%     xlabel('Set of Possible Delays \tau (xT_s)');
%     ylabel('Deviance of Transmitted Pulse');
% end
% 
% %% Delay Plot AWGN - 16 Antenas Scheme
% figure;
% for ii = 1:satAmount;
%     subplot(2,2,ii);
%     plot(LOS_delays/Ts,real(deviance_LOS_est(:,:,ii)));
%     hold on;
%     for snr = [1 6 11 16];
%         plot(LOS_delays/Ts,real(deviance_LOS_S_est_MA_16(:,:,ii,snr)));
%     end
%     hold off;
%     grid on;
%     legend('Noise Free Environment', '-30 dB SNR',...
%             '-20 dB SNR', '-10 dB SNR', '0 dB SNR');
%     xlabel('Set of Possible Delays \tau (xT_s)');
%     ylabel('Deviance of Transmitted Pulse');
% end
% %% DOA Plot
% 
% % figure;
% % subplot(2,2,1);
% % bar(phaseTx(1,:)*180/pi, [1 1 1],0.1);
% % axis([-100 100 0 1.1]);
% % legend('Sem Ru�do');
% % xlabel('�ngulo (Graus)');
% % 
% % subplot(2,2,2);
% % bar(phaseTx(1,:)*180/pi, [1 1 1],0.1);
% % hold on;
% % bar(phaseTx_est(1,:)*180/pi, [1 1 1],0.4,'r');
% % hold off;
% % axis([-100 100 0 1.1]); 
% % legend('Sem Ru�do', '-30 dB');
% % xlabel('�ngulo (Graus)');
% % 
% % subplot(2,2,3);
% % bar(phaseTx(1,:)*180/pi, [1 1 1],0.1);
% % hold on;
% % bar(phaseTx_est(6,:)*180/pi, [1 1 1],0.1,'r');
% % hold off;
% % axis([-100 100 0 1.1]); 
% % legend('Sem Ru�do', '-20 dB');
% % xlabel('�ngulo (Graus)');
% % 
% % subplot(2,2,4);
% % bar(phaseTx(1,:)*180/pi, [1 1 1],0.1);
% % hold on;
% % bar(phaseTx_est(11,:)*180/pi, [1 1 1],0.1,'r');
% % hold off;
% % axis([-100 100 0 1.1]); 
% % legend('Sem Ru�do', '-10 dB');
% % xlabel('�ngulo (Graus)');
% 
% %% DOA Plot 5 antenas
% 
% figure;
% subplot(2,2,1);
% % bar(phaseTx_est(1,:)*180/pi, [1 1 1],0.1);
% % hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_5(1,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% legend('DOA of Tx','4 ULA -30 dB');
% % legend('1 Antena -30 dB', 'DOA de Tx','4 ULA -30 dB');
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,2);
% 
% % bar(phaseTx_est(6,:)*180/pi, [1 1 1],0.1);
% % hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_5(6,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% legend('DOA of Tx','4 ULA -20 dB');
% % legend('1 Antena -20 dB', 'DOA de Tx','4 ULA -20 dB');
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,3);
% 
% % bar(phaseTx_est(11,:)*180/pi, [1 1 1],0.1);
% % hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_5(11,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% legend('DOA of Tx','4 ULA -10 dB');
% %legend('1 Antena -10 dB', 'DOA de Tx','4 ULA -10 dB')
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,4);
% 
% %bar(phaseTx_est(16,:)*180/pi, [1 1 1],0.1);
% %hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_5(16,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% legend('DOA of Tx','4 ULA 0 dB');
% %legend('1 Antena 0 dB','DOA de Tx','4 ULA 0 dB');
% xlabel('Angle (Degrees)');
% 
% %% DOA Plot 16 antenas
% 
% figure;
% subplot(2,2,1);
% %bar(phaseTx_est(1,:)*180/pi, [1 1 1],0.1);
% %hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_16(1,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% legend('DOA of Tx','16 ULA -30 dB');
% %legend('1 Antena -30 dB', 'DOA de Tx','16 ULA -30 dB');
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,2);
% 
% %bar(phaseTx_est(6,:)*180/pi, [1 1 1],0.1);
% %hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_16(6,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% %legend('1 Antena -20 dB', 'DOA de Tx','16 ULA -20 dB');
% legend('DOA of Tx','16 ULA -20 dB');
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,3);
% 
% %bar(phaseTx_est(11,:)*180/pi, [1 1 1],0.1);
% %hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_16(11,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% %legend('1 Antena -10 dB', 'DOA de Tx','16 ULA -10 dB');
% legend('DOA of Tx','16 ULA -10 dB');
% xlabel('Angle (Degrees)');
% 
% subplot(2,2,4);
% 
% %bar(phaseTx_est(16,:)*180/pi, [1 1 1],0.1);
% %hold on;
% bar(phaseTx(1,:)*180/pi, [1 1 1],0.1, 'g');
% hold on;
% bar(phase_MA_est_16(16,:)*180/pi, [1 1 1],0.1,'r');
% hold off;
% axis([-180 180 0 1.1]); 
% %legend('1 Antena 0 dB', 'DOA de Tx','16 ULA 0 dB');
% legend('DOA of Tx','16 ULA 0 dB');
% xlabel('Angle (Degrees)');
% 
% %% RMSE 1 antenna
% figure;
% for ii = 1:4
% subplot(2,2,ii);
% if ii ~= 4
% semilogy(SNR, RMSE(ii,:),'o-.', 'LineWidth', 2, ...
% 'MarkerFaceColor', [0 0 0]);
% hold on;
% else
% subplot(2,2,ii);
% semilogy(SNR, sum(RMSE)/3,'o-', 'LineWidth', 2, ...
% 'MarkerFaceColor', [0 0 0]);
% hold off;
% grid on;
% end
% xlabel('Signal to Noise Ratio (SNR).');
% ylabel('RMSE - Delay Error.');
% legend('No Scheme');
% end
% %% DOA 1 antenna
% figure;
% for ii = 1:4
% subplot(2,2,ii);
% if ii ~= 4
% semilogy(SNR, RMSE_DOA(:,ii),'o-.', 'LineWidth', 2, ...
% 'MarkerFaceColor', [0 0 0]);
% hold on;
% else
% subplot(2,2,ii);
% semilogy(SNR, sum(RMSE_DOA')/3,'o-', 'LineWidth', 2, ...
% 'MarkerFaceColor', [0 0 0]);
% hold off;
% grid on;
% end
% xlabel('Signal to Noise Ratio (SNR).');
% ylabel('RMSE - DOA Error.');
% legend('No Scheme');
% end