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

% calculate the pressure using greens function in free field
% distances from sourc to receivers
R = zeros(M,M);
for i = 1:M
    for ii = 1:M
        R(i,ii) = norm([X(i,ii) Y(i,ii) Z(i,ii)] - xyz);
    end
end

% calculating pressure using greens function in free space
p_ta = 1i*w*rho*Q*exp(1i*(w-k*R))./(4*pi*R);

% Plot
figure
surf(X,Y,20*log10(abs(p_ta)/(20*10^-6)))
colorbar;
shading flat;
zlabel('Sound Pressure Level [dB]');
xlabel('x')
ylabel('y')
title('Task A - SPL of the monopole at (0,0,-0.1)')

% The sound pressure level looks as expected with a high level in the
% center and a decreasing level as the distance to the location of the
% monopole is increased.

figure
contourf(X,Y,angle(p_ta))
shading flat;
hc=colorbar;
title(hc,'Phase [rad]');
xlabel('x')
ylabel('y')
title('Task A - Phase of the pressure at (0,0,-0.1)')

% The phase of the pressure also looks as expected, as it is possible to
% see the phase approximately matches the wavelength. 

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
    surf(KX,KY,abs(squeeze(P_tb(i,:,:))))
    colorbar;
    shading flat;
    zlabel('Pressure [Pa]');
    xlabel('k_x')
    ylabel('k_y')
    title(['Task B - Pressure, f = ' num2str(f(i)) ' Hz at (0,0,-0.1)'])
end

% We decided to plot for 2kHz and 5kHz as well. The plots show very well
% that for the lower frequencies at 500Hz a very simple result is obtained
% as opposed to 2kHz and 5kHz where the wavelength is much shorter. As a
% lower frequency results in a lower wavenumber, k, high frequencies will
% have a large number of propagating waves in z. Lower frequencies, as 500,
% will instead have evanescent waves as k^2 < kx^2+ky^2. Also, there appears 
% aliasing (or folding) at high frequencies because the distance between
% the microphones is not sufficiently small.

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
surf(KX,KY,abs(P_tc))
colorbar;
shading flat;
zlabel('Pressure [Pa]');
xlabel('k_x')
ylabel('k_y')
title('Task C - Pressure, f = 1000 Hz at (L/2,L/2,-0.1)')

% As the monopole is now placed near a corner, it is possible to see the
% shape of the wavenumber spectrum change accordingly. The maximum values
% are seen to move in a bend. A large amount of propagating waves are also
% seen at this point as low kx and ky values are required for evanescent
% waves.

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
surf(KX,KY,abs(P_td))
colorbar;
shading flat;
zlabel('Pressure [Pa]');
xlabel('k_x')
ylabel('k_y')
title('Task D - Pressure, f = 1000 Hz at (0,0,L+3)')

figure
contourf(X,Y,abs(p_td))
shading flat;
hc=colorbar;
title(hc,'Pressure [Pa]');
xlabel('x')
ylabel('y')
title('Task D - Pressure, f = 1000 Hz at (0,0,L+3)')

% The sound pressure is now almost consistent having moved the monopole 3m
% away from the array. This is expected as the array has moved into a far
% field position and the waves are now almost plane resulting in a very
% homogenous distribution of sound pressure. The resulting wave number
% spectrum now consists primarily of propagating waves in z as
% k^2 >= kx^2+ky^2 almost entirely.


%% Task E
% Finally, plot and discuss the wavenumber spectrum of the normal component of the particle velocity,
% when the point source is placed at (x, y, z)= (0, 0 ,-0.1) m.

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
surf(KX,KY,abs(W_te))
colorbar;
shading flat;
zlabel('Particle velocity [m/s]');
xlabel('k_x')
ylabel('k_y')
title('Task E - Particle velocity, f = 1000 Hz at (0,0,-0.1)')

% The normal component of the particle velocity can easily be calculated in
% z_m=0 (in the z origin plane) using Eq (9.218) in the Fundamentals book.
% The behavior resembles the one a monopole should have so our equations
% corretly represent the monopole both in pressure and particle velocity.
