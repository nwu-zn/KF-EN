% Kalman filter of sleep signals sample entropy
clear all
close all

% Data Sampling Rate: 100 Hz
fs=100;

% Read sleep signals
[hdr, record] = edfread('original\SC4001E0.edf');
EEG1=record(1,30630*fs+1:31170*fs);

% [hdr, record] = edfread('original\SC4011E0.edf');
% EEG1=record(1,21540*fs+1:22200*fs);

% [hdr, record] = edfread('original\SC4021E0.edf');
% EEG1=record(1,21870*fs+1:22680*fs);

% [hdr, record] = edfread('original\SC4032E0.edf');
% EEG1=record(1,27630*fs+1:28180*fs);

% Visualization of raw signals
EEG1=EEG1';
lenofEEG1 = length(EEG1);
tEEG = (1:length(EEG1))/fs;
subplot(311)
plot(tEEG, EEG1);
title('EEG');
ylabel('Amplitude (au)')

% Calculate sample entropy
jj=1;
while 1
    window=fs*5;
    if jj*window>lenofEEG1
        break
    end
    sig = EEG1((jj-1)*window+1:jj*window);
    r = 0.15 * std(sig);
    [entropy(jj),envar(jj)] = SampEnVar( 2, r, sig, 1 );
    jj=jj+1;
end

samEn_var=var(entropy);

% Visualization of sample entropy
ttime=(1:length(entropy))*5;
subplot(312)
plot(ttime,entropy)
ylabel('Sample entropy')

% Kalman filter process
A = 1; 
H = 1; 
Q = 0.1;
R = 0.5; 
z = entropy;
x = z(1); 
P  = 1; 
estimated_state = zeros(1, length(entropy));
for k=1:length(entropy) 
    x = A * x;
    P = A * P * A' + Q;
    K = P * H' / (H * P * H' + R);  
    x = x + K * (z(k) - H * x);
    P = (1 - K * H) * P;
	estimated_state(k) = x;
end

% Calculate VRR
estimated_var=var(estimated_state);
VRR = (samEn_var-estimated_var)/samEn_var;

subplot(313)
plot(ttime(1:length(entropy)),entropy)
ylabel('Estimated entropy')

xlabel('Time (s)')
subplot(313)
hold on
plot(ttime(1:length(estimated_state)),estimated_state,'--r')
ylabel('Estimated entropy')
xlabel('Time (s)')
legend('Sample entropy','kalman filter estimation')

