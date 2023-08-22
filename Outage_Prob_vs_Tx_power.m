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

% Define the outage probability
P1=zeros(1,length(pt));
P2=zeros(1,length(pt));

% Define the target rates bps/Hz 
rate1=1;rate2=2;
for u=1:length(pt)
% Estimate the SINR
gamma_1 = a1*pt(u)*g1./(a2*pt(u)*g1+no);
gamma_12 = a1*pt(u)*g2./(a2*pt(u)*g2+no);
gamma_2 = a2*pt(u)*g2/no;


% Estimate the Rates
R1  =  (log2(1+gamma_1));
R12= (log2(1+gamma_12)); 
R2  = (log2(1+gamma_2));

% check for outage
for k=1:N
if R1(k) < rate1    
    P1(u)=P1(u)+1;
end

if R2(k)<rate2 || R12(k)<rate1
P2(u)=P2(u)+1;
end

end
end
PP1 = P1/N;
PP2 = P2/N;

figure
semilogy(Pt_dBm , PP1,'-ok','linewidth',1.5); hold on
semilogy(Pt_dBm , PP2,'--k','linewidth',1.5); hold on
xlabel('Transmit power [dBm]')
ylabel('Outage probability')
grid
legend('User 1 (far user)','User 2 (near user)')
