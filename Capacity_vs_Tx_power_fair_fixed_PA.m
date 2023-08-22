% Initialization
clc; clear all; close all

% Define the variation range of the tranmitted power in dBm
Pt_dBm = 0:30;

% The corresponding value of the transmitted power in linear scale
pt = 1e-3 * db2pow(Pt_dBm);

% Define the number of channel realizzations
N = 1e5;

% Define the distances between the users and base station
d1 = 5000; d2=1000;

% define the pass loss exponent
eta = 4;

% Generate the Rayleigh fading coefficnet for each user with pass loss
h1 = sqrt(d1^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);

% Estimate the channel gain for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

% Define the bandwidth
BW = 1e6;

% Estimate the noise power in watt
no = 1.38064852e-23 * 300 * BW;

% Define the fixed PA coefficnets
b1=0.75; b2=0.25;

% Target Rate bps/Hz
R1 = 1;

% Target SINR
xi = (2^(R1))-1;

for u=1:length(pt)
a1 = xi*(pt(u)*g1+no)./(pt(u)*g1*(xi+1));    
a1(a1>1)=0;
a2=1-a1;

% Sum Rate of fair PA
C1 = log2 (1 + a1.*pt(u).*g1./(a2.*pt(u).*g1+no));
C2 = log2 (1 + a2.*pt(u).*g2/no);

C_sum_fair_PA(u) = mean(C1+C2);


% Sum Rate of fixed PA
C1_ = log2 (1 + b1*pt(u)*g1./(b2*pt(u)*g1+no));
C2_ = log2 (1 + b2*pt(u)*g2/no);

C_sum_fixed_PA(u) = mean(C1_+C2_);     
end

figure
plot(Pt_dBm , C_sum_fair_PA,'-k','linewidth',2); hold on
plot(Pt_dBm , C_sum_fixed_PA,'--k','linewidth',2); hold on
legend('Fair PA','Fixed PA, \alpha_1=0.75, \alpha=0.25')
xlabel('Transmit power [dBm]')
ylabel('Sum rate (bps/Hz)')
grid













