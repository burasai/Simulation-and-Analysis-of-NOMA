%Initialization
clc; clear all; close all

% Define the variation range of the transmitted power in dBm
Pt_dBm = 0 : 40;

% The corresponding transmitted power in linear scale
pt = 10.^((Pt_dBm-30)/10);
%pt = 1e-3 * db2pow(Pt_dBm);

% Define the number of channel realization
N = 1e4;

% Define the distances between the users and the base station
d1 = 1000; d2 = 500;

% Define the pass loss exponent
eta = 4;

% Define the Rayleigh fading channel coefficients for each user with the
% pass loss exponent
h1 = sqrt(d1^-eta)*(randn(N,1)+j*randn(N,1))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(N,1)+j*randn(N,1))/sqrt(2);

% Estimate the channel gain for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

% Define the value of the bandwidth
BW = 1e6;

% Estimate the noise power
no = 1.38064852e-23 * 300 * BW;

% Fixed power allocation coefficients
a1 = 0.75; a2 = 0.25;

% Define the value of the imperfect SIC
err = [0 1e-4 1e-3 1e-2];

for n=1:length(err)
for k=1:length(pt)
R_2_NOMA_d(k) = mean(log2(1 + a2*pt(k)*g2./(err(n)*a1*pt(k)*g2+no)));
end
R_2_NOMA_d_err(n,:) = R_2_NOMA_d;
end

figure
plot(Pt_dBm, R_2_NOMA_d_err(1,:),'-k','linewidth',2); hold on
plot(Pt_dBm, R_2_NOMA_d_err(2,:),'--k','linewidth',2); hold on
plot(Pt_dBm, R_2_NOMA_d_err(3,:),'.-k','linewidth',2); hold on
plot(Pt_dBm, R_2_NOMA_d_err(4,:),':k','linewidth',2); hold on
xlabel('Transmitted power [dBm]')
ylabel('Achievable rate (bps/Hz)')
title('The achievable data rate at user 2 with imperfect SIC')
legend('\epsilon = 0 (perfect SIC)','\epsilon = 10^{-4}','\epsilon = 10^{-3}','\epsilon = 10^{-2}')
grid






