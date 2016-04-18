clear;close all;clc;
load prob1.txt;
load prob2.txt;
for num = 1:250
    freqfea(num) = prob1(3+3*(num-1)+1);
    pin_real(num) = prob1(3+3*(num-1)+2);
    pin_imag(num) = prob1(3+3*(num-1)+3);
    pout_real(num) = prob2(3+3*(num-1)+2);
    pout_imag(num) = prob2(3+3*(num-1)+3);
    
end
pin_spl = 20*log10(abs(pin_real+1i*pin_imag)/(2*10^-5));
pout_spl = 20*log10(abs(pout_real+1i*pout_imag)/(2*10^-5));
% save freqfea.mat freqfea;
% save pin_spl.mat pin_spl;
% save pout_spl.mat pout_spl;

% k = freq*2*pi/340;
subplot(2,1,1)
plot(freqfea, pin_spl);
subplot(2,1,2)
plot(freqfea, pout_spl);