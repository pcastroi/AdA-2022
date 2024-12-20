%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Green's function simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc,
close all
%-------------------------------------------------
% Constants and parameters
%-------------------------------------------------
c = 343; % speed of sound

% room dimension
lx = 4.4; % width
ly = 3.28;  % length
lz = 3.28; % height
V   = lx*ly*lz; % volume

% source and receiver positions
%Position 4
Ps_all=[3.847 0.43 1.288];   %source position
Pr_all=[3.932 3.052 1.366];    %receiver position

%Position 2
% Ps_all=[3.847 0.430 1.288];   %source position
% Pr_all=[2.149 1.981 1.366];   %receiver position

x0 = Ps_all(1,1);% Point of source
y0 = Ps_all(1,2);
z0 = Ps_all(1,3);


x = Pr_all(1,1);
y = Pr_all(1,2);
z = Pr_all(1,3);
% point of interest



% reverberation time 
T_60 = 3.5;

% time constant
tau = T_60/(6*log(10));

% Green function frequencies
f = 50:0.00625:450; 
k = 2*pi*f/c;       % frequency dependent constant k
k_2 = k.^2;         % squaring k
N_samples = length(f);

% ----------------------------------------------------------------------
% Calculate Greens function
% ----------------------------------------------------------------------
% number of modes in each direction
Nx = 24;
Ny = 24;
Nz = 24;
Nx1=Nx;Ny1=Ny;Nz1=Nz;
pos1 = [x,y,z];
Calculate_Greens_Function; % return 'Green'

%-------------------------------------------------
% PLOT simulated Greens function
%-------------------------------------------------
%
G_sim_dB = 20*log10(abs(Green));

% hold on
% plot(f,G_sim_dB,'r');
% hold on
% axis tight
% xlabel('Frequency [Hz]')
% ylabel('Amplitude [dB]')
% title({['Source Position=(',num2str(Ps_all(1)),',',num2str(Ps_all(2)),',',num2str(Ps_all(3)),')'],['Receiver Position=(',num2str(Pr_all(1)),',',num2str(Pr_all(2)),',',num2str(Pr_all(3)),')']})
% grid on


%% Figure parameters

widthofline = 1; % thickness of plot lines
sizeofmarker = 8; % thickness of plot markers (x o ect.)
sizeoffont = 16; % font size for text on plots

% %% Position 1
% 
% %import data - frequency response function 
% F1 = importdata('pos1_Frequency Response H1(p_B_source,p_A_source) - Input.txt');
% F2 = importdata('pos1_Frequency Response H1(Pressure C,p_A_source) - Input.txt');
% F3 = importdata('pos1_Frequency Response H1(Pressure C,p_B_source) - Input.txt');
% 
% %measured frequencies
% ff_me = F1(:,2);
% omega_me = 2*pi*ff_me;
% 
% %measured frequency response function 
% H_re1 = F1(:,3); %real part
% H_im1 = F1(:,4); %imaginary part
% 
% H_re2 = F2(:,3); %real part
% H_im2 = F2(:,4); %imaginary part
% 
% H_re3 = F3(:,3); %real part
% H_im3 = F3(:,4); %imaginary part
% 
% %combining the real and imaginary pary 
% H_mag1 = complex(H_re1,H_im1);
% H_mag2 = complex(H_re2,H_im2);
% H_mag3 = complex(H_re3,H_im3);
% 
% 
% 
% % plot position 1 - source C and A
% figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
% plot(ff_me,20*log10(abs(H_mag2)),f,G_sim_dB,'-', 'LineWidth', widthofline, 'DisplayName', 'this line')
% hold on
% %title('Position 1') % title above plot
% xlabel('f [Hz]') % label for x axis
% ylabel('Green´s function [dB]') % label for y axis
% ylim([-110 10]) % set boundery for y axis
% %xlim([50 250]) % set boundery for x axis
% legend('Measured','Simulated');
% legend('Location','northeast') % place legend on plot
% set(gca,'FontSize',sizeoffont) % change font size
% grid on
% hold off
% 
% 
% %% Position 2
% 
% %import data - frequency response function - position 2
% F12 = importdata('pos2_Frequency Response H1(p_B_source,p_A_source) - Input.txt');
% F22 = importdata('pos2_Frequency Response H1(Pressure C,p_A_source) - Input.txt');
% F32 = importdata('pos2_Frequency Response H1(Pressure C,p_B_source) - Input.txt');
% 
% %measured frequency response function 
% H_re12 = F12(:,3); %real part
% H_im12 = F12(:,4); %imaginary part
% 
% H_re22 = F22(:,3); %real part
% H_im22 = F22(:,4); %imaginary part
% 
% H_re32 = F32(:,3); %real part
% H_im32 = F32(:,4); %imaginary part
% 
% %combining the real and imaginary pary 
% H_mag12 = complex(H_re12,H_im12);
% H_mag22 = complex(H_re22,H_im22);
% H_mag32 = complex(H_re32,H_im32);
% 
% 
% % plot position 2 - source C and A
% figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
% plot(ff_me,20*log10(abs(H_mag22)),f,G_sim_dB,'-', 'LineWidth', widthofline, 'DisplayName', 'this line')
% hold on
% %title('Position 2') % title above plot
% xlabel('f [Hz]') % label for x axis
% ylabel('Green´s function [dB]') % label for y axis
% ylim([-110 10]) % set boundery for y axis
% %xlim([50 250]) % set boundery for x axis
% legend('Measured','Simulated');
% legend('Location','northeast') % place legend on plot
% set(gca,'FontSize',sizeoffont) % change font size
% grid on
% hold off

%% position 3


%import data - frequency response function - position 2
[band13, f13, F13] = read_pulse_2021('Frequency Response H1(p_B_source,p_A_source) - Input.txt');
[band23, f23, F23] = read_pulse_2021('Frequency Response H1(Pressure C,p_A_source) - Input.txt');
[band33, f33, F33] = read_pulse_2021('Frequency Response H1(Pressure C,p_B_source) - Input.txt');

% %measured frequency response function 
% H_re13 = F13(:,3); %real part
% H_im13 = F13(:,4); %imaginary part
% 
% H_re23 = F23(:,3); %real part
% H_im23 = F23(:,4); %imaginary part
% 
% H_re33 = F33(:,3); %real part
% H_im33 = F33(:,4); %imaginary part
% 
% %combining the real and imaginary pary 
% H_mag13 = complex(H_re13,H_im13);
% H_mag23 = complex(H_re23,H_im23);
% H_mag33 = complex(H_re33,H_im33);


% plot position 2 - source C and A
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
plot(f33,20*log10(abs(F33)),f,G_sim_dB,'-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Position 3') % title above plot
xlabel('f [Hz]') % label for x axis
ylabel('Green´s function [dB]') % label for y axis
ylim([-110 10]) % set boundery for y axis
%xlim([50 250]) % set boundery for x axis
legend('Measured','Simulated');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold of


