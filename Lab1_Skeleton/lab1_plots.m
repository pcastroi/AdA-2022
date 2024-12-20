%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021
clear; clc; close all;
% COMPLETE THE FUNCTIONS WITH YOUR CODE, THEN USE THEM TO RUN EXPERIMENTS IN THIS SCRIPT

% DEFINE PARAMETERS FOR YOUR DUCT EXPERIMENT
a = 0.20; % in m
b = 0.50; % in m
x = 0;% in m
y = 0; % in m
z = 0; % in m
x0 = a; % in m
y0 = b; % in m

c = 343;
f = 0:1500;

n_max = 2;
m_max = 5;

[G,cutoff_f] = greens_function_duct(f, x, y, z, m_max, n_max, x0, y0, a, b);
G_dB = 20*log10(abs(G));

plot(f,G_dB);
for i=1:length(cutoff_f(:,1))
    for j=1:length(cutoff_f(1,:))
        if cutoff_f(i,j) <= max(f)
            xline(cutoff_f(i,j),':',['Mode (',num2str(i-1),',',num2str(j-1),')'])
        end
    end
end
hold on

xlabel('Frequency [Hz]')
ylabel('Sound pressure level [dB]')
title('Cut-off frequencies of the Green´s function')
grid on
