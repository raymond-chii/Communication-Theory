clc
clear 
close all

%%%% Lei(Raymond) Chi ps 06

%% c
rolloff = 0.3;
span = 4;
sps = 8; 
Rb = 1e6; 


Rs = Rb / log2(span);

% Design the square-root raised cosine filter using rcosdesign
filterCoeff = rcosdesign(rolloff, span, sps, 'sqrt');
%% d
p = filterCoeff;
g = conv(p, fliplr(p)); 
time_p = 0:length(p)-1;
time_g = 0:length(g)-1;

subplot(2,1,1);
stem(time_p, p, 'Marker', 'o', 'LineWidth', 1.5);
title('Impulse Response p[n]');
xlabel('n');
ylabel('Amplitude');
grid on;
xlim([0, 70]);
subplot(2,1,2);
stem(time_g, g, 'Marker', 'o', 'LineWidth', 1.5);
title('Matched Filter Output g[n]');
xlabel('n');
ylabel('Amplitude');
grid on;
% xlim([0, 35]);

%% e

[maxValue, maxIndex] = max(g)

%% f
T_s = 1 / Rs;
A_0 = 1; 
A_max = max(abs(A_0));
worstISI = sum(abs(g))*A_max; 
sig_Power = A_max^2;
worstISI_SIR_dB = 10 * log10(sig_Power / (worstISI^2))


%% g
W = (1+rolloff)*Rs/2; 

t = linspace(0, span / Rs, span * sps + 1);
f = linspace(0, W, length(t));
F_s = Rs * sps;

G_f = abs(freqz(filterCoeff, 1, f, F_s)).^2;
% g_t = ifft(G_f);

% Plot
figure;
stem(f, G_f, 'LineWidth', 1.5);
title('Impulse Response within [0, W]');
xlabel('frequency');
ylabel('Mag');
grid on;


%% h
SNR_values = [inf, -5, 5, 10];
ava_ber = zeros(length(SNR_values), 10);
ava_ms = zeros(length(SNR_values), 10);

numofBits = 10^6;

for t = 1:10
    for SNR_dB_index = 1:length(SNR_values) 
        SNR_dB = SNR_values(SNR_dB_index);
        
        bits = randi([0, 1], 1, numofBits);
        
        qpskSymbols = zeros(1, numofBits / 2);
        
        for n = 1:length(bits)/2
            p = bits(2*n);
            imp = bits(2*n-1);
            if (imp == 0) && (p == 0)
                qpskSymbols(n) = exp(1j * pi/4);  % 45 degrees
            elseif (imp == 1) && (p == 0)
                qpskSymbols(n) = exp(1j * 3*pi/4);  % 135 degrees
            elseif (imp == 1) && (p == 1)
                qpskSymbols(n) = exp(1j * 5*pi/4);  % 225 degrees
            elseif (imp == 0) && (p == 1)
                qpskSymbols(n) = exp(1j * 7*pi/4);  % 315 degrees
            end
        end
        
        sigma = sqrt(10^(-SNR_dB/10));
        noise = sigma * (randn(size(qpskSymbols)) + 1j * randn(size(qpskSymbols)));
        receivedSignal = qpskSymbols + noise;
        
        filterCoeff = rcosdesign(rolloff, span, sps, 'sqrt');
        filteredSignal = conv(receivedSignal, filterCoeff);
        receivedSymbols = filteredSignal(1:sps:end);
       
        decodedBits = real(receivedSymbols) > 0;
        
        
        bits = bits(1:length(decodedBits));
        numErrors = sum(decodedBits ~= bits);
        ava_ber(SNR_dB_index,t) = numErrors / numofBits;
        
        qpskSymbols = qpskSymbols(1:length(receivedSymbols));
        ava_ms(SNR_dB_index,t) = mean((qpskSymbols - receivedSymbols).^2);

    end
end

disp(['Average BER: ' num2str(mean(ava_ber(:)))]);
disp(['Average RMS: ' num2str(mean(ava_ms(:)))]);

