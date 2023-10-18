clc
clear
close all

%%%% Lei(Raymond) Chi Comp theory ps01

%question 1

cumulative_sum = 0;
E_s_values= zeros(1,4);
k = [2:2:8];
index_forEs = 1;
dmin = 1; 
Eb_plot = zeros(1,4); 



for idx = k

    M = 2.^idx; 
    constellation = zeros(1, M);
    for index = 1:M
        real_part = 2 * (index - 1) - (M - 1);
        imaginary_part = 0;  
        constellation(index) = real_part + 1i * imaginary_part;
    end
        sv = constellation.*constellation; 
        dp = sum(sv);
        var = dp^2;
        cumulative_sum = cumulative_sum + var;
        E_s = cumulative_sum/(2.^idx);
        E_b_values(index_forEs) = E_s
        index_forEs = index_forEs + 1; 
end
figure; 
plot(k, E_b_values, '-o'); 
hold on;
xlabel('k');
ylabel('Eb/dmin^2');
title('Eb/dmin^2 vs. k');
grid on;




hold on;
for idx = 1:numel(k)
    text(k(idx), E_b_values(idx), sprintf('k=%d', k(idx)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
hold off;

%part c 
n = zeros(1, 4);
for idx = 1:4
    M = 2^k(idx); 
    
    Tb = 1 / dmin; 
    Rs = 1 / Tb;
    Rb = k(idx)/Tb; 
    BPS = Rb/Rs;
    BPD = log2(M) / k(idx);
    n(idx) = BPD/BPS;
end

figure; 
plot(k, n, '-o');
hold on;
xlabel('k');
ylabel('n');
title('n vs. k');
grid on;

hold on; 
for idx = 1:4
    text(k(idx), n(idx), sprintf('k=%d', k(idx)), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
hold off; 
%%question 2
T = 0.25; 
%part c
A = 1/T;


%part d 
W = 1; 

f = 0:0.01:5;


Y = exp((-log(2)/2)*(f/W).^2-(1j*pi*f*T)).*sinc(f*T)*A*T;
magY = abs(Y);

magY_dB = 20*log10(magY); 

figure;
plot(f, magY_dB);
hold on;
axis([0 5 -60 0])
xlabel('Frequency (f)');
ylabel('|Y(f)| (dB)');
title('Y(f) vs. f');
grid on;

%part e

for i = 1:length(magY_dB)
    if magY_dB(i) <= -50
        B0 = f(i);  
        break; 
    end
end
B0


%part f

t = -0.4:0.01:0.8;
y = A.* qfunc((2*pi*W * (t - T)) / sqrt(log(2)));
figure;
plot(t, y);
hold on;
xlabel('time');
ylabel('output y(t)');
title('y(t) vs time');
grid on;


%part f-5
T0 = 0;
for i = 1:length(y)
    if y(i) <= 4*0.1
        T0 = t(i);  
        break; 
    end
end
T0

