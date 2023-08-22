% Initialization
clc; clear all; close all

% Define the transmitted power in dBm
Pt_dBm = 30;

% The corresponding value of the transmitted power in linear scale
pt = 1e-3 * db2pow(Pt_dBm);

% Define the number of channel realizations
N = 1e5;

% Define the distances between the users and base station
d1 = 1000; d2 = 500;

% Define the pass loss exponent
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

% Fixed PA coefficients
a1 = 0.75; a2 = 0.25;

% Define the far user target rate (R*)
r = 0.5 : 0.5 : 10;

% Define the outage probabiltiy
far_user_outage_prob_FPA = zeros(1,length(r));
near_user_outage_prob_FPA = zeros(1,length(r));
far_user_outage_prob_DPA = zeros(1,length(r));
near_user_outage_prob_DPA = zeros(1,length(r));

for u=1:length(r)
% Define the SINR of the far user    
xi=(2^(r(u)))-1;
aa1 = min(1,(xi*(pt*g1+no)./(pt*g1*(xi+1))));
aa2=1-aa1;

% Capacity of fixed PA
C_f = log2 (1 + a1*pt*g1./(a2*pt*g1+no));
C_nf = log2 (1 + a1*pt*g2./(a2*pt*g2+no));
C_n = log2 (1 + a2*pt*g2./(no));

% Capacity of dynamic PA
C_f_ = log2 (1 + aa1.*pt.*g1./(aa2.*pt.*g1+no));
C_nf_ = log2 (1 + aa1.*pt.*g2./(aa2.*pt.*g2+no));
C_n_ = log2 (1 + aa2.*pt.*g2/(no));

% Estimation of the outage Prob. of fixed PA
for k=1:N
if C_f(k)<r(u)
    far_user_outage_prob_FPA(u) = far_user_outage_prob_FPA(u)+1;
end

if C_n(k) < r(u)  ||  C_nf(k) < r(u)
    near_user_outage_prob_FPA(u) = near_user_outage_prob_FPA(u) + 1;
end

% Estimation of the outage Prob. of dynamic PA
if C_f_(k)<r(u)
    far_user_outage_prob_DPA(u) = far_user_outage_prob_DPA(u)+1;
end

if aa1(k) == 0

if C_n_(k) < r(u)
    near_user_outage_prob_DPA(u) = near_user_outage_prob_DPA(u) + 1;
end      

else

if C_n_(k) < r(u)  ||  C_nf_(k) < r(u)
    near_user_outage_prob_DPA(u) = near_user_outage_prob_DPA(u) + 1;
end    

end
    
    
end
end

P1 = far_user_outage_prob_FPA/N;
P2 = near_user_outage_prob_FPA/N;
P3 = far_user_outage_prob_DPA/N;
P4 =near_user_outage_prob_DPA/N;

figure
plot(r , P1 , '--+r','linewidth',2); hold on
plot(r , P2 , '--ob','linewidth',2); hold on
plot(r , P3 , 'r','linewidth',2); hold on
plot(r , P4 , 'b','linewidth',2); hold on
xlabel('Target rate of far user (R*) bps/Hz')
ylabel('Outage probability')
title('P = 30 dBm')
legend('Far user (fixed PA)','Near user (fixed PA)','Far user (fair PA)','Near user (fair PA)')
grid
axis([0.5 10 0 1])














