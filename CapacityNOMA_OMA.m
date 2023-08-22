% Initialization
clc; clear all; close all;

% Define the variation range of the transmitted power in dB
Pt_dB = -114 : -54;

% The coressponding transmitted power in linear scale
pt = 10.^(Pt_dB/10);
%pt = db2pow(Pt_dB);

% Define the number of channel realization
N = 1e4;

% Define the distances between the users and base station
d1=6; d2=3; d3=2;

% Define the pass loss exponent
eta = 4;

%Generate the Rayleigh fading channel coefficients for each user with pass loss exponent
h1 = sqrt(d1^-eta)*(randn(N,1) + j*randn(N,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N,1) + j*randn(N,1))/sqrt(2);
h3 = sqrt(d3^-eta)*(randn(N,1) + j*randn(N,1))/sqrt(2);

% Estime the channel gains for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;

% Defien the value of the bandwidth
BW=1e9;

% Estime the noise power
no = 1.38064852e-23 * 300 * BW; % in watt
No = 10*log10(no); % in dB

SNR_dB = Pt_dB - No;

% First Fixed power allocation coefficients
%a1 = 0.75; a2 = a1 * (1-a1); a3 = 1-(a1+a2);

% Second Fixed power allocation coefficients
d= (d1^2)+(d2^2)+(d3^2);
a1 = (d1^2)/d; a2 = (d2^2)/d; a3 = (d3^2)/d;

% Inintialization of the capacity calculation
C_NOMA = zeros(1, length(pt));
C_OMA = zeros(1, length(pt));

for k=1:length(pt)
% Downlink NOMA Capacity Calculation
C_NOMA_1_d = log2(1 + a1*pt(k)*g1./(a2*pt(k)*g1 + a3*pt(k)*g1 + no));% User 1
C_NOMA_2_d = log2(1 + a2*pt(k)*g2./(a3*pt(k)*g2+no)); % User 2
C_NOMA_3_d = log2(1 + a3*pt(k)*g3/no); % User 3
C_NOMA_sum_d(k) = mean(C_NOMA_1_d+C_NOMA_2_d+C_NOMA_3_d);

% Downlink/Uplink OMA Capacity Calculation
C_OMA_1_u_d = (1/3)*log2(1 + pt(k)*g1/no); % User 1
C_OMA_2_u_d = (1/3)*log2(1 + pt(k)*g2/no); % User 2
C_OMA_3_u_d = (1/3)*log2(1 + pt(k)*g3/no); % User 3
C_OMA_sum_u_d(k) = mean(C_OMA_1_u_d + C_OMA_2_u_d + C_OMA_3_u_d);

% Uplink NOMA Capacity Calculation
C_NOMA_1_u = log2(1 + a1*pt(k)*g1./(a2*pt(k)*g2+a3*pt(k)*g3+no));%User 1
C_NOMA_2_u = log2(1 + a2*pt(k)*g2./(a3*pt(k)*g3+no));%User 2
C_NOMA_3_u = log2(1 + a3*pt(k)*g3/no);
C_NOMA_sum_u(k) = mean(C_NOMA_1_u+C_NOMA_2_u+C_NOMA_3_u);
end
figure
plot(SNR_dB , C_NOMA_sum_d, '-k', 'linewidth',2); hold on
plot(SNR_dB , C_OMA_sum_u_d, '--k', 'linewidth',2); hold on
xlabel('SNR [dB]');
ylabel('Achievable rate (bps/Hz)');
legend('NOMA','OMA')
title('Capacity of downlink NOMA vs. OMA')
grid
axis([0 60 0 16])

figure
plot(SNR_dB , C_NOMA_sum_u, '-r', 'linewidth',2); hold on
plot(SNR_dB , C_OMA_sum_u_d, '--r', 'linewidth',2); hold on
xlabel('SNR [dB]');
ylabel('Achievable rate (bps/Hz)');
legend('NOMA','OMA')
title('Capacity of uplink NOMA vs. OMA')
grid
axis([0 60 0 16])

