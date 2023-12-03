clc
clear
close all

%%%% Lei(Raymond) Chi ps05 


%% question 1
% a
P_e = 10^(-5);
p = -0.217;
gemma_b1 = (qfuncinv(P_e))^2 / (1 - real(p));
gemma_b1_db = 10 * log10(gemma_b1) 
% b
gemma_b2 = (qfuncinv(P_e))^2;
gemma_b2_db = 10 * log10(gemma_b2)
% c 
gemma_b3 = -log(2*P_e)*2;
gemma_b3_db = 10 * log10(gemma_b3)


%% question 2

% "00" = 45 deg
% "01" = 135 deg
% "11" = -135 deg
% "10" = -45 deg

first = 0 + 45;
second = first - 135;
third = second - 135;
forth = third + 45; 
fifth = forth + -45; 

encode_phase = [first, second, third, forth, fifth]

%% question 3

% c

E_0 = 1;
d_1 = 2*sqrt(E_0) % d_min
d_2 = sqrt(2)*2*sqrt(E_0)
d_3 = sqrt(4*sqrt(E_0)^2+(2*sqrt(E_0))^2)
d_4 = sqrt((-sqrt(E_0))^2+(-sqrt(E_0))^2)
d_5 = 3*2*sqrt(E_0)
d_6 = 2*2*sqrt(E_0)


distances = [d_1, d_2, d_3, d_4, d_5, d_6];
sorted_distances = sort(distances);

counters = 0;

for i = 1:length(sorted_distances)
    for j = 1:length(distances)
        if distances(j) == sorted_distances(i)
            counters = counters + 1;
        end
    end
end

disp(['Total Pairs: ' num2str(counters)]);
% d
M = 8;
gamma_b_values = linspace(0, 10, 100);  
beta = 0.5:0.5:4;
a_values = ones(1, M); 

P_b_bound = zeros(size(gamma_b_values));
for i = 1:length(gamma_b_values)
    P_b_bound(i) = sum(a_values .* qfunc(sqrt(beta .* gamma_b_values(i))));
end

% e
max_value = max(beta);
min_value = min(beta);
P_b_approx_max_at_2dB = a_values(1) * qfunc(sqrt(max_value * 2));
P_b_approx_max_at_8dB = a_values(1) * qfunc(sqrt(max_value * 8));

P_b_approx_min_at_2dB = a_values(1) * qfunc(sqrt(min_value * 2));
P_b_approx_min_at_8dB = a_values(1) * qfunc(sqrt(min_value * 8));

% f
gamma_b_values_dB = 2:0.1:8;
gamma_b_values = 10.^(gamma_b_values_dB / 10);

P_err_bound = zeros(length(gamma_b_values), M);
for i = 1:length(gamma_b_values)
    P_err_bound(i, :) = a_values .* qfunc(sqrt(beta .* gamma_b_values(i)));
end
P_err_bound_sum = sum(P_err_bound, 2);


ratio_at_2dB = (1 - P_b_approx_max_at_2dB) ./ P_err_bound_sum(gamma_b_values_dB == 2);
ratio_at_8dB = (1 - P_b_approx_max_at_8dB) ./ P_err_bound_sum(gamma_b_values_dB == 8);

disp(['Ratio at 2dB: ' num2str(ratio_at_2dB)]);
disp(['Ratio at 8dB: ' num2str(ratio_at_8dB)]);

figure;
hold on;
semilogy(gamma_b_values_dB, 1 - P_err_bound_sum, 'DisplayName', 'Union Bound');
semilogy(gamma_b_values_dB, ones(size(gamma_b_values_dB)) .* (1 - P_b_approx_max_at_8dB), '-.', 'DisplayName', 'Approximation (Max at 8dB)');
semilogy(gamma_b_values_dB, ones(size(gamma_b_values_dB)) .* (1 - P_b_approx_min_at_2dB), ':', 'DisplayName', 'Approximation (Min at 2dB)');
semilogy(gamma_b_values_dB, ones(size(gamma_b_values_dB)) .* (1 - P_b_approx_min_at_8dB), '-+', 'DisplayName', 'Approximation (Min at 8dB)');
hold off;
xlabel('Bit SNR (gamma_b) [dB]');
ylabel('1 - Probability of Bit Error (1 - P_b)');
title('Union Bound and Approximations for 8-ary QAM with Gray Coding');
legend('Location', 'Best');
grid on;
