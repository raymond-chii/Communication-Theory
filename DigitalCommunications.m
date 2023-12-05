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

d1 = sqrt(4*E_0);
d2 = sqrt(8*E_0);
d3 = sqrt(16*E_0);
d4 = sqrt(20*E_0);
d5 = sqrt(36*E_0);
d6 = sqrt(40*E_0);

points = [-3, 1; -1, 1; 1, 1; 3, 1; -3, -1; -1, -1; 1, -1; 3, -1];


% distances = [d1, d2, d3, d4, d5, d6];
% sorted_distances = sort(distances);

count_d1 = 0;  
count_d2 = 0; 
count_d3 = 0;  
count_d4 = 0; 
count_d5 = 0;  
count_d6 = 0; 

distance = 0; 

for i = 1:size(points, 1)
    for j = 1:size(points, 1)
        distance = norm(points(i, :) - points(j, :));
        if distance == d1
            count_d1 = count_d1 + 1;
        elseif distance == d2
            count_d2 = count_d2 + 1;
        elseif distance == d3
            count_d3 = count_d3 + 1;
        elseif distance == d4
            count_d4 = count_d4 + 1;
        elseif distance == d5
            count_d5 = count_d5 + 1;
        elseif distance == d6
            count_d6 = count_d6 + 1;
        end
    end
end
counters = count_d1 + count_d2 + count_d3 + count_d4 + count_d5 + count_d6;

disp(['Total Pairs: ' num2str(counters)]);
% d
% M = 8;
gamma_b_values_dB = 2:0.1:8;
gamma_b_values = 10.^(gamma_b_values_dB / 10);

beta = [2/2, 4/2, 8/2, 10/2, 18/2, 20/2];
a_values = [count_d1, count_d2, count_d3, count_d4, count_d5, count_d6]/24;

P_b_bound = zeros(size(gamma_b_values_dB));
for i = 1:length(a_values)
    P_b_bound = P_b_bound + (a_values(i) * qfunc(sqrt(beta(i) * gamma_b_values)))/3;
end

% max_value = max(beta);
% min_value = min(beta);
% P_b_approx_max_at_2dB = a_values(1) * qfunc(sqrt(max_value * 2));
% P_b_approx_max_at_8dB = a_values(1) * qfunc(sqrt(max_value * 8));
% 
% P_b_approx_min_at_2dB = a_values(1) * qfunc(sqrt(min_value * 2));
% P_b_approx_min_at_8dB = a_values(1) * qfunc(sqrt(min_value * 8));]

% f
ratioat2dB = 1 - (P_b_bound(1)/sum(P_b_bound))
ratioat8dB = 1 - (P_b_bound(end)/sum(P_b_bound))
figure;
semilogy(gamma_b_values_dB, P_b_bound, 'DisplayName', 'Union Bound');
xlabel('Bit SNR / gamma_b (dB)');
ylabel('P_b');
title('Union Bound and Approximations for 8-ary QAM with Gray Coding');
legend('Location', 'Best');
grid on;