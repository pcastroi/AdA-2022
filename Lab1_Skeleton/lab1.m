%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021
% clear; clc; close all;
% COMPLETE THE FUNCTIONS WITH YOUR CODE, THEN USE THEM TO RUN EXPERIMENTS IN THIS SCRIPT

% DEFINE PARAMETERS FOR YOUR DUCT EXPERIMENT
 a = 0.2; % in m
 b = 0.25; % in m
 x = 0.2;% in m
 y = 0.2; % in m
 z = 0.1; % in m
 z2 = 1;
 x0 = 0; % in m
 y0 = 0; % in m
 n_max = 20;
 m_max = 20;
 
 c = 343;
 f = 0:5000;
 

G = greens_function_duct(f, x, y, z, m_max, n_max, x0, y0, a, b);
G_2 = greens_function_duct(f, x, y, z2, m_max, n_max, x0, y0, a, b);
% figure
G_sim_dB = 20*log10(abs(G));
G_sim_dB2 = 20*log10(abs(G_2));
plot(f,G_sim_dB);
hold on
plot(f,G_sim_dB2,'--r');
xlabel('Frequency [Hz]')
ylabel('Sound pressure level [dB]')

z_2 = linspace(0,1,100);
f_1 = 800;
f_2 = 3200;
G_3 = zeros(length(z_2));
G_4 = zeros(length(z_2));
for ii = 1:length(z_2)
    G_3(ii) = greens_function_duct(f_1, x, y, z_2(ii), m_max, n_max, x0, y0, a, b);
    G_4(ii) = greens_function_duct(f_2, x, y, z_2(ii), m_max, n_max, x0, y0, a, b);
end
figure
G_sim_dB3 = 20*log10(abs(G_3));
G_sim_dB4 = 20*log10(abs(G_4));
plot(z_2,G_sim_dB3);
hold on
plot(z_2,G_sim_dB4,'--r');
xlabel('Distance [m]')
ylabel('Sound pressure level [dB]')

