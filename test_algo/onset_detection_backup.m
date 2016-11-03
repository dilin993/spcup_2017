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

s = log(MF*abs(s)); % log of mel spectrum of the signal

[M, N] = size(s);

ds = zeros(size(s));
ds(:,2:N) =  diff(s,1,2)*fs;

ds(ds<0) = 0;

acc = sum(ds,1);

acc_mean = smooth(acc,3);
acc_mean = acc_mean';
plot(t,acc);

hold on;

plot(t,acc_mean);

n_cur = acc - acc_mean;
n_cur(n_cur<0)= 0;
TH = mean(n_cur)+ 4*std(n_cur);
n_cur(n_cur<TH) = 0;
hold on;

plot(t,n_cur);


figure;

t_all = linspace(0,(L-1)/fs,L);
plot(t_all, y);

hold on;


t_mark = t(n_cur>0)';
marks = zeros(size(t_mark,1),1);
stem(t_mark,marks);

player.play()
