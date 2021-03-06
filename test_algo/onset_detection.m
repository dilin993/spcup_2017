close all;
clear all;

[y, fs] = audioread('..\training_set\open\open_001.wav');

L = size(y, 1);
WINDOW = 1024;
OVERLAP = 512;
NFFT = 2^(nextpow2(WINDOW)+2);
MF = melfb(40,NFFT,fs);
player = audioplayer(y, fs);

[s,f,t] = spectrogram(y,WINDOW,OVERLAP,NFFT,fs);

s = 10*log(MF*abs(s)); % mel spectrum of the signal in dB

N = size(s,2);

ds = zeros(size(s));
ds(:,2:N) =  diff(s,1,2)*fs;

ds(ds<0) = 0;

acc = sum(ds,1);

figure;
plot(t,acc);
title('Raw onset signal');
xlabel('time(s)')

p = dcblock(0.4,fs);
b = [1 -1];
a = [1 -p];
acc0 = filtfilt(b,a,acc);
u = mean(acc0);
sigma = std(acc0);
acc0 = (acc0-u)/sigma;

figure;
plot(t,acc0);
title('normalized onset signal');
xlabel('time(s)')

sigma = 5;
size = 800;
x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

O = conv(acc0,gaussFilter,'same');

figure;
plot(t,O);
title('smoothed onset signal');
xlabel('time(s)')



