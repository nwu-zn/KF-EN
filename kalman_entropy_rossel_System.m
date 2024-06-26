% kalman filter of rossel system signals sample entropy
clear all
close all

% Data Sampling Rate: 100 Hz
fs = 100;

% Generate rossel signals
a=0.38;
b=0.2;
c=5.7;		
x(1)=1;
y(1)=0;
z(1)=0;
h=0.01;
for i=1:50000
    K1=-(y(i)+z(i));
    L1=x(i)+a*y(i);
    M1=b+x(i)*z(i)-c*z(i);

    K2=-((y(i)+h/2*L1)+(z(i)+h/2*M1));
    L2=(x(i)+h/2*K1)+a*(y(i)+h/2*L1);
    M2=b+(x(i)+h/2*K1)*(z(i)+h/2*M1)-c*(z(i)+h/2*M1);

    K3=-((y(i)+h/2*L2)+(z(i)+h/2*M2));
    L3=(x(i)+h/2*K2)+a*(y(i)+h/2*L2);
    M3=b+(x(i)+h/2*K2)*(z(i)+h/2*M2)-c*(z(i)+h/2*M2); 

    K4=-((y(i)+h*L3)+(z(i)+h*M3));
    L4=(x(i)+h*K3)+a*(y(i)+h*L3);
    M4=b+(x(i)+h*K3)*(z(i)+h*M3)-c*(z(i)+h*M3);  

    x(i+1)=x(i)+h/6*(K1+2*K2+2*K3+K4);
    y(i+1)=y(i)+h/6*(L1+2*L2+2*L3+L4);
    z(i+1)=z(i)+h/6*(M1+2*M2+2*M3+M4);
end
grid;
x=x';
y=y';
z=z';

% Visualization of raw signals
EEG1 = x;
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




