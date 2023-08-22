% Initialization
clc; clear all ;

% Define the tranmsitting power in dBm
Pt_dBm = 30;

% The corresponding value in linear scale
pt = 1e-3 * db2pow(Pt_dBm);

% Difine the value of the bandwidth
BW = 1e4;

% Estimate the noise power in watt
no = 1.38064852e-23 * 300 * BW;

% Define the pass loss exponent
eta = 4;

% Define the number of user
K = 2 : 2 : 20;

%  Defien the number of channel realization
N = 1e3;

% Define the number of Monte Carlo simulation
mc = 1e3;


alpha1 = 0.6 : 0.05 : 0.95;

for c=1:length(alpha1)

for e=1:mc

for u = 1 :length (K)
% Define the distances for each user
dist = randi ([ 200 1500], K(u), 1);

% Sort users in descending order based on distance
d = sort(dist, 'descend');

% Generate The Rayleigh fading coefficient of zeor mean and unit variance
% for all users
H = (randn(K(u),N)+j*randn(K(u),N))/sqrt(2);
    
% Introducing the pass loss exponent
h = bsxfun(@times, sqrt(d.^-eta), H);
    
% Estimate the channel gains
g = (abs(h)).^2;
    
% Power Allocation
PA = zeros (K(u), 1);

for m=1:K(u)
    if m == K(u)
PA(m) = 1-sum(PA(1:m-1));
    else
        PA(m) = (1- sum(PA(1:m-1)))*alpha1(c);
    end
end

% The sum term in the denominator
PA_sum = zeros (K(u),1);
for m=1:K(u)
PA_sum(m) = sum(PA(m+1:end));
end

% Estimating the achievable rate capacity
C = log2 (1+ bsxfun(@times, PA,g)*pt./(bsxfun(@times, PA_sum, g)*pt+no));
C_sum (e, u) = sum(mean(C.'));
end
end

C_sum_NOMA(c,:) = mean(C_sum);
end

%figure
plot(K, C_sum_NOMA(1,:),'-k','linewidth',2);hold on
plot(K, C_sum_NOMA(2,:),'--k','linewidth',2);hold on
plot(K, C_sum_NOMA(3,:),':k','linewidth',2);hold on
plot(K, C_sum_NOMA(4,:),'-r','linewidth',2);hold on
plot(K, C_sum_NOMA(5,:),'--r','linewidth',2);hold on
plot(K, C_sum_NOMA(6,:),':r','linewidth',2);hold on
plot(K, C_sum_NOMA(7,:),'-b','linewidth',2);hold on
plot(K, C_sum_NOMA(8,:),'--b','linewidth',2);hold on
xlabel('Number of users')
ylabel('Sum rate (bps/Hz)')
title('Pt = 30 dBm')
legend('\alpha_1=0.6','\alpha_1=0.65','\alpha_1=0.7','\alpha_1=0.75','\alpha_1=0.8','\alpha_1=0.85','\alpha_1=0.9','\alpha_1=0.95')
grid













