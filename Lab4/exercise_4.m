%% Lab 4 - sound power determination
addpath 'lab 4 - measurements'\
%% Figure parameters

widthofline = 2; % thickness of plot lines
sizeofmarker = 8; % thickness of plot markers (x o ect.)
sizeoffont = 16; % font size for text on plots
savepictures = 0; % if 1 saves plots as epsc. if 0 dont save

%% Constant parameter

c = 343; %speed of sound [m/s]
rho = 1.2; %densoty of air
V = 245; %room volume [m^3]
S = 240; %surface area of the room [m^2]

P0 = 1.0*10^-12; %reference sound power value

correction1 = [-0.1 -0.2 -0.6 -1.0 -1.6]; %dB
correction1_freq = [4000 5000 63000 8000 10000];

correction2 = [1.4 1.4 1.4 1.7 1.8 1.8 2.2 2.4 3 3.3 3.6]; %dB
correction2_freq = [1000 1250 1600 2000 2500 3150 4000 5000 63000 8000 10000];



%% Mean reverberation time

%import data - mean reverberation time
T_mea = importdata('Mean Reverberation Time Spectrum.txt');

%measured frequencies
f_mea = T_mea(:,2);

%measured reverberation time
T60 = T_mea(:,3);

% %plot
% figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
% semilogx(f_mea,T60, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
% hold on
% title('Mean reverberation time') % title above plot
% xlabel('f [Hz]') % label for x axis
% xticks([f_mea])
% ylabel('T [s]') % label for y axis
% %ylim([-100 -10]) % set boundery for y axis
% %legend('');
% %legend('Location','northeast') % place legend on plot
% set(gca,'FontSize',sizeoffont) % change font size
% grid on
% hold off


figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
bar(log10(f_mea),T60)
hold on
%title('Mean reverberation time') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f_mea])
ylabel('T [s]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
%legend('');
%legend('Location','northeast') % place legend on plot
%xticks([63 125 250 500 1000 2000 4000 8000])
set(gca,'xtick',log10([63 125 250 500 1000 2000 4000 8000]));
set(gca,'xticklabel',10.^get(gca,'xtick'));
%set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off


%% Measured p_rms^2 - mean square pressure - dipole

%close to wall/boundary

%import data 
p1_wall = importdata('Dipole_Wall_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f1_mea = p1_wall(:,2);

%measured mean square sound pressure
p1 = p1_wall(:,3);

% %plot
% figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
% semilogx(f1_mea,p1, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
% hold on
% %title('Mean square sound pressure') % title above plot
% xlabel('f [Hz]') % label for x axis
% %xticks([f1_mea])
% ylabel('Mean square pressure [Pa]') % label for y axis
% %ylim([-100 -10]) % set boundery for y axis
% %legend('');
% %legend('Location','northeast') % place legend on plot
% set(gca,'FontSize',sizeoffont) % change font size
% grid on
% hold off


%out in the room

%import data - mean reverberation time
p11_far = importdata('Dipole_Room1_Autospectrum(Signal 1) - Input.txt');
p12_far = importdata('Dipole_Room2_Autospectrum(Signal 1) - Input.txt');
p13_far = importdata('Dipole_Room3_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f11_mea = p11_far(:,2);

%measured mean square sound pressure
p11 = p11_far(:,3);
p12 = p12_far(:,3);
p13 = p13_far(:,3);

%mean of the three source positions
p1_out = mean([p11 p12 p13],2);


% %plot
% figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
% semilogx(f11_mea,p1_out, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
% hold on
% %title('Mean square sound pressure') % title above plot
% xlabel('f [Hz]') % label for x axis
% %xticks([f1_mea])
% ylabel('Mean square pressure [Pa]') % label for y axis
% %ylim([-100 -10]) % set boundery for y axis
% %legend('');
% %legend('Location','northeast') % place legend on plot
% set(gca,'FontSize',sizeoffont) % change font size
% grid on
% hold off


%combined plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f1_mea,20*log10(p1),f11_mea,20*log10(p1_out), '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('Mean square pressure [dB]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the walls','Far from the walls');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off



%% Sound power, dipole

lambda = c./f1_mea; %wavelength

P0 = 1.0*10^-12; %reference sound power value

%far from the wall
P_dipole = 13.8*((V.*p1_out)./(rho*c^2*T60)).*(1+(S.*lambda/8*V));

P_dipole_dB = 10*log10(P_dipole/P0);

%close to the wall/boundary
P_dipole_wall = 13.8*((V.*p1)./(rho*c^2*T60)).*(1+(S.*lambda/8*V));

P_dipole_wall_dB = 10*log10(P_dipole_wall/P0);

%plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f11_mea,P_dipole_dB,f1_mea,P_dipole_wall_dB, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('sound power level [dB]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the wall','Far from the wall');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off





%% Measured p_rms^2 - mean square pressure - low impedance source

%close to wall/boundary

%import data 
p2_wall = importdata('Low_Wall_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f2_mea = p2_wall(:,2);

%measured mean square sound pressure
p2 = p2_wall(:,3);


%out in the room

%import data - mean reverberation time
p21_far = importdata('Low_Room1_Autospectrum(Signal 1) - Input.txt');
p22_far = importdata('Low_Room2_Autospectrum(Signal 1) - Input.txt');
p23_far = importdata('Low_Room3_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f21_mea = p21_far(:,2);

%measured mean square sound pressure
p21 = p21_far(:,3);
p22 = p22_far(:,3);
p23 = p23_far(:,3);

%mean of the three source positions
p2_out = mean([p21 p22 p23],2);


%combined plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f2_mea,p2,f21_mea,p2_out, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('Mean square pressure [Pa]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the wall','Far from the wall');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off


%% Sound power, Low impedance source

lambda1 = c./f2_mea; %wavelength

%far from the wall
P_low = 13.8*((V.*p2_out)./(rho*c^2*T60)).*(1+(S.*lambda1/8*V));

P_low_dB = 10*log10(P_low/P0);

%close to the wall/boundary
P_low_wall = 13.8*((V.*p2)./(rho*c^2*T60)).*(1+(S.*lambda1/8*V));

P_low_wall_dB = 10*log10(P_low_wall/P0);

%plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f21_mea,P_low_dB,f2_mea,P_low_wall_dB, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('Acoustical sound power [W]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the wall','far from the wall');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off


%% Measured p_rms^2 - mean square pressure - High impedance source

%close to wall/boundary

%import data 
p3_wall = importdata('High_Wall_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f3_mea = p3_wall(:,2);

%measured mean square sound pressure
p3 = p3_wall(:,3);


%out in the room

%import data - mean reverberation time
p31_far = importdata('High_Room1_Autospectrum(Signal 1) - Input.txt');
p32_far = importdata('High_Room2_Autospectrum(Signal 1) - Input.txt');
p33_far = importdata('High_Room3_Autospectrum(Signal 1) - Input.txt');

%measured frequencies
f31_mea = p31_far(:,2);

%measured mean square sound pressure
p31 = p31_far(:,3);
p32 = p32_far(:,3);
p33 = p33_far(:,3);

%mean of the three source positions
p3_out = mean([p31 p32 p33],2);



%combined plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f3_mea,p3,f31_mea,p3_out, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('Mean square pressure [Pa]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the wall','Far from the wall');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off


%% Sound power, High impedance source

lambda1 = c./f3_mea; %wavelength

%far from the wall
P_high = 13.8*((V.*p3_out)./(rho*c^2*T60)).*(1+(S.*lambda1/8*V));

P_high_dB = 10*log10(P_high/P0);

%close to the wall/boundary
P_high_wall = 13.8*((V.*p3)./(rho*c^2*T60)).*(1+(S.*lambda1/8*V));

P_high_wall_dB = 10*log10(P_high_wall/P0);

%plot
figure('Name','Name of figure','NumberTitle','off','pos',[100 100 1200 700]);
semilogx(f31_mea,P_high_dB,f3_mea,P_high_wall_dB, '-', 'LineWidth', widthofline, 'DisplayName', 'this line')
hold on
%title('Mean square sound pressure') % title above plot
xlabel('f [Hz]') % label for x axis
%xticks([f1_mea])
ylabel('Acoustical sound power [W]') % label for y axis
%ylim([-100 -10]) % set boundery for y axis
legend('Close to the wall','far from the wall');
legend('Location','northeast') % place legend on plot
set(gca,'FontSize',sizeoffont) % change font size
grid on
hold off



