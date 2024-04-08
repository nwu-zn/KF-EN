%Kalman filter of logistic map signals entropy
clear all
close all

% Data Sampling Rate: 100 Hz
fs = 100;

% Generate logistic map signals
r = 3.9; 
x0 = 0.5; 
n = 1000; 
noise_length = 50000; 
x = zeros(1, n + noise_length);
x(1) = x0; 
for i = 2:(n + noise_length)
    x(i) = r * x(i-1) * (1 - x(i-1));
end
noise = x((n+1):end);

% Visualization of raw signals
EEG1 = noise;
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
estimated_var=var(estimated_state);
VRR = (samEn_var-estimated_var)/samEn_var;


subplot(313)
plot(ttime(1:length(entropy)),entropy)
ylabel('Estimated sample entropy')

xlabel('Time (s)')
subplot(313)
hold on
plot(ttime(1:length(estimated_state)),estimated_state,'--r')
ylabel('Estimated sample entropy')
xlabel('Time (s)')
legend('Sample entropy','kalman filter estimation')
