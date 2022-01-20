
clear; close all; clc;


Fs = 10e3;
t = 0:1/Fs:3;

y = 2.*sin(2*pi*30*t) + 2.*sin(2*pi*150*t);

wind = 200;


[Pxx,f] = pwelch(y,[],0,[],Fs,'power','onesided');


plot(f(1:100),pow2db(Pxx(1:100)))

avp1 = mean(Pxx);

N = length(y);
avp2 = bandpower(y,Fs,(N-1)*[0,Fs/(2*N)]);
    
pRMs = rms(y)^2;


avp3 = (1/N) * sum(abs(y).^2);