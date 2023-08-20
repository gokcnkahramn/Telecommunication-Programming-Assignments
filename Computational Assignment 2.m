%% Question 1
% For each question part, you need to run codes "Run Section"
% First of all, we enter our parameter values ​​in this section.
valuePam = [4, 8, 16]; % They are different PAM levels
listSNR = [5:0.5:13, 8:0.5:17, 10:0.5:22]; % SNR list in the unit of dB
% Then, I need to convert SNR to linear
convertSNR = 10.^(listSNR./10);
% Numeric bit error of transmissions
transNum = 100000;
% I arranged P errors on both theoretical and simulated.
theoPerror = zeros(length(valuePam), length(convertSNR));
sim_Perror= zeros(length(valuePam), length(convertSNR));

for m = 1:length(valuePam)
    numPam = valuePam(m); % This is total numbitER of PAM levels    
    for s = 1:length(convertSNR)
        SNR_last = convertSNR(s); % This is last SNR which is converted        
        % This is our random symbols 
        rand_sym = randi([0, numPam-1], 1, transNum);        
        % It makes PAM demodulation
        mod_pam = 2*rand_sym - (numPam-1);      
        % This part cretes AWGN
        noise = sqrt(1/(2*SNR_last)) * randn(1, transNum);
        sym_rec = mod_pam + noise;        
        % It makes PAM demodulation
        decSym = (sym_rec + (numPam-1)) / 2;    
        % This function counts symbol error
        sym_err = sum(decSym ~= rand_sym);        
        % This is simulated probability of symbol error
        sim_Perror(m, s) = sym_err / transNum;        
        % This is theoretical probability of symbol error
        theoPerror(m, s) = 2 * (1 - 1/numPam) * qfunc(sqrt(3 * log2(numPam) * SNR_last / (numPam^2 - 1)));
    end
end

% This part makes plotting
figure;
semilogy(listSNR, sim_Perror(1,:), 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(listSNR, theoPerror(1,:), 'r-', 'LineWidth', 2);
semilogy(listSNR, sim_Perror(2,:), 'gd-', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(listSNR, theoPerror(2,:), 'm-', 'LineWidth', 2);
semilogy(listSNR, sim_Perror(3,:), 'ks-', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(listSNR, theoPerror(3,:), 'c-', 'LineWidth', 2);
hold off;
grid on;
xlabel('SNR/bit in  (dB)');
ylabel('Probability of Error');
legend('M = 4 Simulated', 'M = 4 Theoretical', 'M = 8 Simulated', 'M = 8 Theoretical', 'M = 16 Simulated', 'M = 16 Theoretical');
title('Probability of Error for each PAM');

%% Question 2a

% This is parameters part
bit_num = 1000;
T = 1;
Tsampling = T/20; % Sampling period
listSNR = -10:1:20;
% Then, I need to convert SNR to linear
convertSNR = 10.^(listSNR./10);
% Select pulse shape (example: raised cosine pulse)
t = -5*T:Tsampling:5*T;
roll_off_factor = 0.5; % Roll-off factor
pulse = (sin(pi*t/T) ./ (pi*t/T)) .* (cos(roll_off_factor*pi*t/T) ./ (1 - (2*roll_off_factor*t/T).^2));
% Matched filter
filt_match = fliplr(pulse);
% Initialize error rate arrays
theoPerror = zeros(1, length(convertSNR));
sim_Perror= zeros(1, length(convertSNR));

for s = 1:length(convertSNR)
    SNR_last = convertSNR(s); % converted SNR    
    % This is generated random bits
    bits = randi([0, 1], 1, bit_num);
    % This function makes PAM modulation
    mod_pam = 2*bits - 1;
    % This makes Upsample PAM symbols
    sym_upsamp = zeros(1, length(mod_pam) * (T/Tsampling));
    sym_upsamp(1:T/Tsampling:end) = mod_pam; 
    % We need to transmit those signal
    sig_trans = conv(sym_upsamp, pulse, 'same');
    % We need to AWGN for this case
    stdNoise = 1 / sqrt(2*SNR_last);
    noise = stdNoise * randn(1, length(sig_trans));
    sig_rec = sig_trans + noise;
    % Output for matched sampled filter
    out_match = conv(sig_rec, filt_match, 'same');
    % Output for sampled filter
    out_samp = out_match(1:T/Tsampling:end);   
    % This function makes PAM demodulation
    decSym = (out_samp > 0);
    % Bit error count
    err_bit = sum(decSym ~= bits);
    % Probability of error in simulated
    sim_Perror(s) = err_bit / bit_num; 
    % Probability of error in theoretical
    theoPerror(s) = qfunc(sqrt(SNR_last));
end

% This part is for plotting
figure;
semilogy(listSNR, sim_Perror, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(listSNR, theoPerror, 'r-', 'LineWidth', 2);
hold off;

grid on;
xlabel('SNR in dB');
ylabel('Probability of Error');
legend('Simulated', 'Theoretical');
title('Error Rates for Binary PAM');
%% Question 2b

% Enter parameters
bit_num = 1000; 
T = 1; 
Tsampling = T/20; 
var_jitter = [0.1, 0.3, 0.5]; % List of jitter variances
SNRdB = 10; % SNR in dB
% Then, I need to convert SNR to linear
convertSNR = 10^(SNRdB/10);
% I chose to build pulse shape
t = -5*T:Tsampling:5*T;
roll_off_factor = 0.5;
pulse = (sin(pi*t/T) ./ (pi*t/T)) .* (cos(roll_off_factor*pi*t/T) ./ (1 - (2*roll_off_factor*t/T).^2));
% I need to match my filter
filt_match = fliplr(pulse);
% Then I initialized error on array
sim_Perror = zeros(1, length(var_jitter));

for j = 1:length(var_jitter)
    sum_var_jitter = var_jitter(j); % Sum of variance of jitter filter
    total_bit_errors = 0;
    
    for n = 1:bit_num
        % random bits generated
        bit = randi([0, 1]);
        % This is PAM modulation
        sym_pam = 2*bit - 1;
        
        % This is upsampled symbol
        sym_upsampled = zeros(1, length(pulse));
        sym_upsampled(1:T/Tsampling:end) = sym_pam;
        % This is translated signal
        sig_trans = conv(sym_upsampled, pulse, 'same');    
        % Then, I added AWGN
        stdNoise = 1 / sqrt(2*convertSNR);
        noise = stdNoise * randn(1, length(sig_trans));
        sig_rec = sig_trans + noise; 
        % Jitter added to sampled localization
        jitter = sqrt(sum_var_jitter) * randn(1);
        ind_sampled = round(n * T / Tsampling) + round(jitter / Tsampling);
        % This is our matched filter output with jitter
        out_match = sym_upsampled .* filt_match;
        % This is our sample matched filter output
        out_samp = out_match(max(1, min(length(out_match), ind_sampled)));
        % This is PAM demodulation
        decoded_bit = (out_samp > 0);
        % This part will calculate bit error
        bit_error = xor(bit, decoded_bit);
        total_bit_errors = total_bit_errors + bit_error;
    end    
    % This is simulation of total error
    sim_Perror(j) = total_bit_errors / bit_num;
end

% This part is for plotting
figure;
semilogy(var_jitter, sim_Perror, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Jitter Variance');
ylabel('Probability of Error');
title('Simulated Error Rates with Jitter');
%% Question 3

% This part is for our parameters
numPam = [4, 8, 16];  % Constellation sizes
bps = log2(numPam);
range_snr = -10:0.5:22;
% Theoretical Bit Error Probability for every case
theo_bitER4 =  2 * (1 - 1/sqrt(4)) * (1/2) * (1 - erf(sqrt(10.^(range_snr / 10) * bps(1) / (2*(4-1)))));  % Example for M = 4
theo_bitER8 = 2 * (1 - 1/sqrt(8)) * (1/2) * (1 - erf(sqrt(10.^(range_snr / 10) * bps(2) / (2*(8-1)))));  % Calculate upper bound for M = 8
theo_bitER16 = 2 * (1 - 1/sqrt(16)) * (1/2) * (1 - erf(sqrt(10.^(range_snr / 10) * bps(3) / (2*(16-1)))));  % Calculate upper bound for M = 16
% Symbol error rate and bit error rate fixing
symboER = zeros(length(numPam), length(range_snr));
bitER = zeros(length(numPam), length(range_snr));
% Repeat for each sizes
for m = 1:length(numPam)
    % Calculate the numbitER of bits
    cou_bit = ceil(10000 / bps(1)) * bps(1);
    % Repeat for different SNR values
    for indexSNR = 1:length(range_snr)
        % Firstly, I generated random symbols
        rand_sym = randi([0 numPam(m)-1], 1, floor(cou_bit / bps(m)));
        % Modulation happened
        phase = rand_sym * (2*pi/numPam(m));
        transSig = exp(1j * phase);
        % I added AWGN channel
        SNR = 10^(range_snr(indexSNR) / 10);
        noise_power = 1 / (SNR * bps(m));
        noise = sqrt(noise_power/2) * (randn(size(transSig)) + 1j * randn(size(transSig)));
        recSig = transSig + noise;
        % Demodulation happened
        phase_demod = angle(recSig);
        sym_demod = round(phase_demod * numPam(m) / (2*pi));
        % Symbol and Bit error sum
        err_symbol = sum(sym_demod ~= rand_sym);
        % For bit error calculation, forms fixed
        mat_sym = reshape(rand_sym, bps(m), []);
        demodMat_sym = reshape(sym_demod, bps(m), []);
        error_bit = sum(mat_sym ~= demodMat_sym, 'all');
        symboER(m, indexSNR) = err_symbol / length(rand_sym);
        bitER(m, indexSNR) = error_bit / (length(rand_sym) * bps(m));
    end
end

% Plotting
figure;
semilogy(range_snr, symboER(1,:), 'o-', 'DisplayName', 'symboER (M = 4)');
hold on;
semilogy(range_snr, bitER(1,:), 'o-', 'DisplayName', 'M = 4 Bit Error Rate');
semilogy(range_snr, bitER(2,:), 'o-', 'DisplayName', 'M = 8 Bit Error Rate');
semilogy(range_snr, bitER(3,:), 'o-', 'DisplayName', 'M = 16 Bit Error Rate');
semilogy(range_snr, theo_bitER4, '--', 'DisplayName', 'M = 4 Theoretical Bit Error Rate');
semilogy(range_snr, theo_bitER8, '--', 'DisplayName', 'M = 8 Theoretical Bit Error Rate');
semilogy(range_snr, theo_bitER16, '--', 'DisplayName', 'M = 16 Theoretical Bit Error Rate');
xlabel('SNR per Bit');
ylabel('Probability of Error');
title('Symbol Error Rate and Bit Error Rate Performance');
legend('Location', 'best');
grid on;
%% Question 4

% This is our parameters
numPam = [4, 8, 16]; 
bps = log2(numPam);
range_snr = -10:2:20; 
% Arrange PSK and PAM results
ser_pam = zeros(length(numPam), length(range_snr));
ber_pam = zeros(length(numPam), length(range_snr));
ser_psk = zeros(length(numPam), length(range_snr));
ber_psk = zeros(length(numPam), length(range_snr));
% Try for all M values
for btw = 1:length(numPam)
    % Repeat for all SNR values
    for indexSNR = 1:length(range_snr)
        % Random symbols generated
        rand_sym = randi([0 numPam(btw)-1], 1, 10000);
        % PAM
        mod_pam = 2*rand_sym - (numPam(btw)-1);
        % PSK
        sym_psk = exp(1j * 2*pi*rand_sym/numPam(btw));
        % SNR calculation
        snr = 10^(range_snr(indexSNR)/10);
        noise_var_pam = var(mod_pam) / (2 * snr);
        noise_var_psk = var(sym_psk) / (2 * snr);
        % Gaussian noise generation
        noise_pam = sqrt(noise_var_pam) * randn(size(mod_pam));
        noise_psk = sqrt(noise_var_psk) * randn(size(sym_psk));
        % Received signals with noise
        rec_sym_pam = mod_pam + noise_pam;
        rec_sym_psk = sym_psk + noise_psk;
        % Demodulation for PAM
        dec_sym_pam = (rec_sym_pam + (numPam(btw)-1)) / 2;
        % Demodulation for PSK
        dec_sym_psk = angle(rec_sym_psk) * (numPam(btw) / (2*pi));
        dec_sym_psk = mod(round(dec_sym_psk), numPam(btw));
        % Symbol Error Calculation
        error4pam = sum(dec_sym_pam ~= rand_sym);
        error4psk = sum(dec_sym_psk ~= rand_sym);
        % Bit Error Calculation
        biterror4pam = sum(dec_sym_pam ~= rand_sym) * bps(btw);
        biterror4psk = sum(dec_sym_psk ~= rand_sym) * bps(btw);
        % Calculate Symbol Error Rate and Bit Error Rate
        ser_pam(btw, indexSNR) = error4pam / length(rand_sym);
        ber_pam(btw, indexSNR) = biterror4pam / (length(rand_sym) * bps(btw));
        ser_psk(btw, indexSNR) = error4psk / length(rand_sym);
        ber_psk(btw, indexSNR) = biterror4psk / (length(rand_sym) * bps(btw));
    end
end

% This part helps us to compare Question 1 and 3
fprintf('Comparison of PAM and PSK:\n');
for btw = 1:length(numPam)
    fprintf('\nFor M = %d:\n', numPam(btw));
    for indexSNR = 1:length(range_snr)
        fprintf('SNR: %.2f\n', range_snr(indexSNR));
        fprintf('PAM Symbol Error: %.4e\n', ser_pam(btw, indexSNR));
        fprintf('PAM Bit Error: %.4e\n', ber_pam(btw, indexSNR));
        fprintf('PSK Symbol Error: %.4e\n', ser_psk(btw, indexSNR));
        fprintf('PSK Bit Error: %.4e\n', ber_psk(btw, indexSNR));
        % This part makes which one dominates
        if ber_pam(btw, indexSNR) < ber_psk(btw, indexSNR)
            fprintf('PAM is better than PSK.\n');
        elseif ber_pam(btw, indexSNR) > ber_psk(btw, indexSNR)
            fprintf('PSK is better than PAM.\n');
        else
            fprintf('PAM and PSK have similar performance.\n');
        end
    end
    fprintf('---Loading---\n');
end

% Plotting
figure;
colors = {'b', 'r', 'g'};
markers = {'o', 's', 'd'};
for btw = 1:length(numPam)
    semilogy(range_snr, ser_pam(btw, :), [colors{btw} '-' markers{btw}], 'DisplayName', sprintf('PAM symboER (M = %d)', numPam(btw)));
    hold on;
    semilogy(range_snr, ber_pam(btw, :), [colors{btw} '--' markers{btw}], 'DisplayName', sprintf('PAM bitER (M = %d)', numPam(btw)));
    semilogy(range_snr, ser_psk(btw, :), [colors{btw} ':' markers{btw}], 'DisplayName', sprintf('PSK symboER (M = %d)', numPam(btw)));
    semilogy(range_snr, ber_psk(btw, :), [colors{btw} '-.' markers{btw}], 'DisplayName', sprintf('PSK bitER (M = %d)', numPam(btw)));
end
xlabel('SNR per Bit');
ylabel('Probability of Error');
title('Comparison of PAM and PSK');
legend('Location', 'best');
grid on;
hold off;