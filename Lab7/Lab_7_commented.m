clear; clc; close all;

xyz = [0 0 -0.1]; % location of monopole
Q = 10^-4; % volume velocity of monopole

c = 343; % speed of sound in air
rho = 1.2; % density
f = 1000; % frequency radiated by monopole
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
lambda = c/f; % wavelength

L = 1; % array dimensions
M = 15; % number of microphones in each direction
z_m = 0;

%% Task A
[X,Y,Z] = meshgrid(-L/2:(L/(M-1)):L/2,-L/2:(L/(M-1)):L/2,z_m);

% calculate the pressure using greens function in free field?

% distances from source to receivers
R = zeros(M,M);
for i = 1:M
    for ii = 1:M
        R(i,ii) = norm([X(i,ii) Y(i,ii) Z(i,ii)] - xyz);
    end
end

p_ta = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R); % calculating pressure using greens function in free space

% Plot
figure
contourf(X,Y,20*log10(abs(p_ta)/(20*10^-6)),'Color','k')
shading flat;
colorbar;
xlabel('x')
ylabel('y')
zlabel('Sound Pressure Level [dB]')
grid on
title('Task A')



figure
contourf(X,Y,angle(p_ta),'Color','k')
shading flat;
colorbar;
xlabel('x')
ylabel('y')
zlabel 'Phase [rad]'
grid on
title('Task A')



%% Task B
% Calculate the wavenumber transform of this sound field, and examine the resulting wavenumber
% spectrum. Analyze what happens if the source radiates at a frequency of 500 Hz, instead of 1 kHz.
f=[500 1000 2000 5000];
N=M*10;
kx = linspace(-(2*pi)/(2*L/(M-1)),(2*pi)/(2*L/(M-1)),N);
ky = kx;
[KX,KY] = meshgrid(kx,ky);

for i=1:length(f)
    w = 2*pi*f(i); % angular frequency
    k = w/c; % wavenumber
    kz(i,:) = sqrt(k^2-kx.^2-ky.^2);

    p_tb(i,:,:) = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R);
    
    % Wavenumber transform
    P_tb(i,:,:) = fftshift(fft2(squeeze(p_tb(i,:,:)),N,N));
    
    % Plot
    figure
    contourf(KX,KY,abs(squeeze(P_tb(i,:,:))))
    colorbar;
    xlabel('k_x')
    ylabel('k_y')
    title(['Task B, f = ' num2str(f(i)) ' Hz'])
end
colorbar;
shading flat;

%% Task C
% Examine the wavenumber spectrum of the sound pressure (again at 1 kHz) when the monopole
% source is placed near a corner of the array, for example (x, y, z)= (L/2, L/2 ,-0.1) m.
f=1000;
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
kz = sqrt(k^2-kx.^2-ky.^2);
xyz = [L/2,L/2,-0.1];
% distances from source to receivers
R = zeros(M,M);
for i = 1:M
    for ii = 1:M
        R(i,ii) = norm([X(i,ii) Y(i,ii) Z(i,ii)] - xyz);
    end
end
p_tc = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R);
P_tc = fftshift(fft2(p_tc,N,N));

% Plot
figure
contourf(KX,KY,abs(P_tc))
colorbar;
shading flat;
xlabel('k_x')
ylabel('k_y')
title('Task C')



%% Task D
% Place the monopole 3 m away from the array , and examine both the resulting sound pressure field
% and its wavenumber spectrum
f=1000;
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
kz = sqrt(k^2-kx.^2-ky.^2);
xyz = [0,0,L+3]; % Placing the source 3 meters away in y axis
R = zeros(M,M);
for i = 1:M
    for ii = 1:M
        R(i,ii) = norm([X(i,ii) Y(i,ii) Z(i,ii)] - xyz);
    end
end
p_td = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R);
P_td = fftshift(fft2(p_td,N,N));

% Plot
figure
contourf(KX,KY,abs(P_td))
shading flat;
colorbar;
xlabel('k_x')
ylabel('k_y')
title('Task D')

figure
contourf(X,Y,abs(p_td))
shading flat;
colorbar;
xlabel('x')
ylabel('y')
title('Task D')


%% Task E
% Finally, plot and discuss the wavenumber spectrum of the normal component of the particle velocity,
% when the point source is placed at (x, y, z)= (0, 0 ,-0.1) m. 
% Plot W(kx,ky) -> Eq (9.218) in book.

f=1000;
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
kz = sqrt(k^2-kx.^2-ky.^2);
xyz = [0,0,-0.1]; % Placing the source 3 meters away in y axis
R = zeros(M,M);
for i = 1:M
    for ii = 1:M
        R(i,ii) = norm([X(i,ii) Y(i,ii) Z(i,ii)] - xyz);
    end
end
p_te = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R);
P_te = fftshift(fft2(p_te,N,N));

W_te = P_te.*kz/(k*rho*c);
% Plot
figure
contourf(KX,KY,abs(W_te))
shading flat;
colorbar;
xlabel('k_x')
ylabel('k_y')
title('Task E')