% Initialization
clc ; clear all; close all;

% Define the distancesbetween the users and base station
d1 = 10; d2 = 9; d3 = 4; d4 = 3;

% Define the pass loss exponent
eta = 4;

% Define the number of channel realizations
N = 1e4;

% Generate the Rayleigh fading channel with pass loss exponent
h1 = sqrt(d1^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h3 = sqrt(d3^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h4 = sqrt(d4^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);

% Estime the channel gain for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;
g3 = (abs(h3)).^2;
g4 = (abs(h4)).^2;

% Define the SNR variation in dB
SNR_dB = 20:1:40;

% The corresponding SNR in linear scale
SNR=db2pow(SNR_dB);

% Define the power allocation coefficients
a1=0.75; a2=1-a1;

b1=0.75; b2=b1*(1-b1); b3=b1*(1-(b1+b2)); b4=1-(b1+b2+b3);

for u=1:length(SNR_dB)
% Capacity of N-F paring
C11(u)= mean(log2(1+a1*g1*SNR(u)./(a2*g1*SNR(u)+1)));
C44(u)= mean(log2(1+a2*g4*SNR(u)));
C22(u)= mean(log2(1+a1*g2*SNR(u)./(a2*g2*SNR(u)+1)));
C33(u)= mean(log2(1+a2*g3*SNR(u)));

% Capacity of N-N, F-F paring
C1(u)= mean(log2(1+a1*g1*SNR(u)./(a2*g1*SNR(u)+1)));
C2(u)= mean(log2(1+a2*g2*SNR(u)));
C3(u)= mean(log2(1+a1*g3*SNR(u)./(a2*g3*SNR(u)+1)));
C4(u)= mean(log2(1+a2*g4*SNR(u)));

% Capacity of NOMA
N1(u) = mean(log2(1 + b1*g1*SNR(u)./((b2+b3+b4).*g1*SNR(u)+1)));
N2(u) = mean(log2(1 + b2*g2*SNR(u)./((b3+b4).*g2*SNR(u)+1)));
N3(u) = mean(log2(1 + b3*g3*SNR(u)./((b4).*g2*SNR(u)+1)));
N4(u) = mean(log2(1 + b4*g4*SNR(u)));
    
%Capacity of TDMA
T1(u) = 0.25*mean(log2(1+g1*SNR(u)));    
T2(u) = 0.25*mean(log2(1+g2*SNR(u)));    
T3(u) = 0.25*mean(log2(1+g3*SNR(u)));    
T4(u) = 0.25*mean(log2(1+g4*SNR(u)));          
end

R_sum_N_F = 0.5*(C11+C22+C33+C44);
R_sum_N_N = 0.5*(C1+C2+C3+C4);
R_sum_NOMA=N1+N2+N3+N4;
R_sum_TDMA=T1+T2+T3+T4;

figure
plot(SNR_dB, R_sum_N_F,'-xk','linewidth',1.5); hold on
plot(SNR_dB, R_sum_N_N,'--^k','linewidth',1.5); hold on
plot(SNR_dB, R_sum_NOMA,':*b','linewidth',1.5); hold on
plot(SNR_dB, R_sum_TDMA,'-or','linewidth',1.5); hold on
xlabel('SNR [dB]')
ylabel('Sum rate (bps/Hz)')
legend('Hybrid NOMA N-F paring','Hybrid NOMA N-N, F-F paring','SC-NOMA','TDMA')
grid











