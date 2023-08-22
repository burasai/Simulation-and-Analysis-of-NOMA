% Initialization
clc; clear all; close all;

% Define the number of channel realiazaions
N = 1e5;

% Define the distances between the users and base station
d1 = 1000; d2 = 500;

% Define the fixed power allocation coefficient
a1 = 0.75 ; a2 = 0.25;

% Define the pass loss exponent
eta=4;

% Generate the Rayeligh fading channel for each user with pass loss
h1 = sqrt(d1^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);

% Estimate the gain for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

% Define the transmitting power range in dBm
Pt_dBm = 0 : 2 : 40;

% The corresponding transmitting power in linear scale
pt = 1e-3 * db2pow(Pt_dBm);

% Define the bandwidth
BW = 1e6;

% Estimate the noise power in watt
no = 1.38064852e-23 * 300 * BW;

for u=1:length(pt)
% Estimate the SINR
gamma_1 = a1*pt(u)*g1./(a2*pt(u)*g1+no);
gamma_12 = a1*pt(u)*g2./(a2*pt(u)*g2+no);
gamma_2 = a2*pt(u)*g2/no;


% Average SINR
gamma_1_av(u)  = 10*log10(mean(gamma_1));
gamma_12_av(u)= 10*log10(mean(gamma_12));
gamma_2_av(u)  = 10*log10(mean(gamma_2));

% Estimate the Rates
R_1(u)  =  mean(log2(1+gamma_1));
R_12(u)= mean(log2(1+gamma_12)); 
R_2(u)  = mean(log2(1+gamma_2));
end
figure
plot(Pt_dBm , gamma_1_av,'-k','linewidth',1.5); hold on
plot(Pt_dBm , gamma_12_av,'--k','linewidth',1.5); hold on
plot(Pt_dBm , gamma_2_av,':k','linewidth',1.5); hold on
xlabel('Transmit power [dBm]')
ylabel('Achievable SINR (dB)')
grid
legend('\gamma_1','\gamma_{12}','\gamma_2')

figure
plot(Pt_dBm , R_1,'-r','linewidth',1.5); hold on
plot(Pt_dBm , R_12,'--r','linewidth',1.5); hold on
plot(Pt_dBm , R_2,':r','linewidth',1.5); hold on
xlabel('Transmit power [dBm]')
ylabel('Achievable capacity (dB)')
grid
legend('R_1','R_{12}','R_2')




