clear all; close all; clc;
p=1.2;
c=343;
f=0:10000;
k=(2*pi*f)/343;
Q=10^-3;
h=0.01;
% h=linspace(0.0001,1,length(f));

Pa_dipole=(p*c*h.^2.*k.^4*abs(Q)^2)/(6*pi);
Pa_quadrupole=(2*p*c*h.^4.*k.^6*abs(Q)^2)/(5*pi);

figure
plot(f,Pa_dipole,f,Pa_quadrupole)
set(gca, 'XScale', 'log','YScale','log')
legend('dipole','quadripole')
xlabel('Frequency (Hz)')
ylabel('Sound Power (dB)')
grid on