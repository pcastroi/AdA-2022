%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021
% clear; clc; close all;
% COMPLETE THE FUNCTIONS WITH YOUR CODE, THEN USE THEM TO RUN EXPERIMENTS IN THIS SCRIPT

% DEFINE PARAMETERS FOR YOUR DUCT EXPERIMENT
a = 0.20; % in m
b = 0.20; % in m
x = 0;% in m
y = 0; % in m
z = 0; % in m
x0 = a; % in m
y0 = b; % in m

c = 343;
f = 0:5000;
modes=5;

figure
for ii=0:modes
m_max = 1;
n_max = ii;
G = greens_function_duct(f, x, y, z, m_max, n_max, x0, y0, a, b);
G_dB = 20*log10(abs(G));

plot(f,G_dB);
hold on
end
xlabel('Frequency [Hz]')
ylabel('Sound pressure level [dB]')
title('Green´s function with m_{max} = 1 and n_{max} = [0 1 2 3 4 5]')
legend('(1,0)','(1,1)','(1,2)','(1,3)','(1,4)','(1,5)')
grid on
hold off
