% Initialization
clc; clear all; 

% Define the value of the transmitted power in dBm
Pt_dBm = 0 : 2 : 40; 

% The corresponding value of the tranmitted power in linear scale
pt = 1e-3 * db2pow(Pt_dBm);

% Define the number of channel realization
N = 1e5;

% Define the distances between users and base station
d1 = 1000; d2 = 500;

% Define the pass loss exponent
eta = 4;

% Generate the Rayleigh fading coefficients for each user with the loss
% exponent
h1 = sqrt(d1^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);
h2 = sqrt(d2^-eta)*(randn(1,N)+j*randn(1,N))/sqrt(2);

% Estimate the channel gain for each user
g1 = (abs(h1)).^2;
g2 = (abs(h2)).^2;

% Define the bandwidth
BW = 1e6;

% Estimate the noise power in watt
no = 1.38064852e-23 * 300 * BW;

% Generate the npoise samples for both users
w1 = sqrt(no)*(randn(1,N)+j*randn(1,N))/sqrt(2);
w2 = sqrt(no)*(randn(1,N)+j*randn(1,N))/sqrt(2);

% Generate a uniform distributed data vector of length N
data1 = randi([0 1], 1 , N);
data2 = randi([0 1], 1 , N);

% Mapped to NRZ (1's and -1's)
x1 = 2*data1-1;
x2 = 2*data2-1;

% Define the power allocation coefficient
a1 =0.75; a2 = 1-a1;

for u=1:length(pt)
% assisgned the power allocation to each user and built the NOMA signal
x = sqrt(pt(u))*(sqrt(a1)*x1 + sqrt(a2)*x2);% NOMA signal

% The the AWGN vector and the Rayleigh channel effect
y1 = h1.*x+w1;
y2 = h2.*x+w2;

% Equalization procedure at each user
x1_eq = y1./h1;
x2_eq = y2./h2;

% Direct decoding to recover x1 signal
x1_est = zeros(1, N);
x1_est(x1_eq>0) = 1; % if x1_eq is greater than 0, then set x1_est = 1

% Estimate the BER performance for the first user
BER1(u) = biterr(data1 , x1_est)/N;

% Dorect decoding of x1 at the second user
x12_est = zeros(1, N);
x12_est(x2_eq>0) = 1; % if x1_eq is greater than 0, then set x12_est = 1
x12_est = 2*x12_est - 1;


% x12_est = ones(1, N);
% x12_est(x1_eq<0) = -1; 

y2_est = x2_eq - sqrt(a1) * sqrt(pt(u)) * x12_est;
x2_est = zeros(1, N);
x2_est(y2_est>0) = 1; 

% Estimate the BER performance for the second user
BER2(u) = biterr(data2 , x2_est)/N;
end

figure
semilogy(Pt_dBm , BER1,'-ok','linewidth',1.5); hold on
semilogy(Pt_dBm , BER2,'-^b``','linewidth',1.5); hold on
xlabel('Transmit power [dBm]')
ylabel('BER')
legend('User 1 (far user)','User 2 (near user)')
grid












