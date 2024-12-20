%% Ad. A - Lab 6 - Radiation from a spherical source
clear; clc; close all;

c = 343; % speed of sound
rho = 1.2; % density
n = 150; % Number of samples

theta = linspace(0,2*pi,n);
f = linspace(0,200,n);
k = f*2*pi./c;
w = 2*pi*f;
Am=1;
pnorm_calc = zeros(3,n);
pnorm_calc2 = zeros(3,n);
press = zeros(3,n);

%% Task 1 z=1
z=1;
for m = 0:2
    % spherical
    h = hankelsph(m,z);
    
    % Legendre function of order m
    P = legendre(m,cos(theta));
    
    p_t1(m+1,:)=Am*h.*P(1,:);
end

figure()
polarplot(theta,20*log10(abs(p_t1(1,:))/abs(p_t1(1,1))));
hold on
polarplot(theta,20*log10(abs(p_t1(2,:))/abs(p_t1(2,1))));
polarplot(theta,20*log10(abs(p_t1(3,:))/abs(p_t1(3,1))));
legend
rlim([-30 5])

%% Task 2+3
theta=0;
ka=0.1;
r_a= linspace(1,1000,n*100);
z=ka*r_a;
for m = 0:2
    % spherical
    h = hankelsph(m,z);
    
    % Legendre function of order m
    P = legendre(m,cos(theta));
    
    % derivative Hankel
    dh_t4 = hankelder(m,z);

    p_t2(m+1,:)=h.*P(1,:);
    u_t3(m+1,:)=-1/(1i*rho*c)*Am*P(1,:)*dh_t4;
end

% Normalizing by the far field
figure
plot(r_a,20*log10(abs(p_t2(1,:))/abs(p_t2(1,end)))) 
hold on
plot(r_a,20*log10(abs(p_t2(2,:))/abs(p_t2(2,end))))
plot(r_a,20*log10(abs(p_t2(3,:))/abs(p_t2(3,end))))
set(gca, 'XScale', 'log')
grid on

% Phase angle between particle velocity and pressure
pha_angle=rad2deg(unwrap(angle(u_t3./p_t2)));
figure
plot(r_a,pha_angle)
set(gca, 'XScale', 'log')
grid on
legend

%% Task 4 + 5
nmod=100; % Number of modes
theta = linspace(0,2*pi,1000);
r_a=100;
ka=[0.1 1 5 10];
a=1;
p_ff=zeros([length(ka) length(theta)]);
p_ax=zeros([length(ka) length(theta)]);
p_if=zeros([length(ka) length(theta)]);
for i=1:length(ka)
    z=ka(i)*r_a;
    for m=0:nmod
        h = hankelsph(m,z);
        P = legendre(m,cos(theta));
        dh_t4 = ka(i)*hankelder(m,ka(i)); % This changed
        p_ff(i,:)=p_ff(i,:)+(m+0.5)*h.*P(1,:).*dh_t4.^(-1);
        p_ax(i,:)=p_ax(i,:)+(m+0.5)*h.*dh_t4.^(-1);
    end
    p_if(i,:)=(2*pi*a^2*exp(-1i*ka(i)*r_a))/(4*pi*r_a);
end

p_t4=p_ff./p_ax;
p_t5=p_ff./p_if;

figure
ax = polaraxes;
polarplot(theta,20*log10(abs(p_t4(1,:))));
hold on
polarplot(theta,20*log10(abs(p_t4(2,:))));
polarplot(theta,20*log10(abs(p_t4(3,:))));
polarplot(theta,20*log10(abs(p_t4(4,:))));
legend
ax.ThetaZeroLocation = 'top';
rlim([-40 inf])

figure
ax = polaraxes;
polarplot(theta,20*log10(abs(p_t5(1,:))));
hold on
polarplot(theta,20*log10(abs(p_t5(2,:))));
polarplot(theta,20*log10(abs(p_t5(3,:))));
polarplot(theta,20*log10(abs(p_t5(4,:))));
legend
ax.ThetaZeroLocation = 'top';
rlim([-40 inf])

%% Optional Task
z=1;
theta = linspace(0,2*pi,1000);
phi = linspace(0,2*pi,1000);

for n = 0:2
    for m = 0:2
        % Associated Legendre function of order m,n
        P = legendre(m,cos(theta));
        Y = sqrt((2*m+1)/(4*pi)*factorial(m-n)/factorial(m+n))*P.*exp(1i*n.*phi);
    end
end

figure()
polarplot(theta,20*log10(abs(p_t1(1,:))/abs(p_t1(1,1))));
hold on
polarplot(theta,20*log10(abs(p_t1(2,:))/abs(p_t1(2,1))));
polarplot(theta,20*log10(abs(p_t1(3,:))/abs(p_t1(3,1))));
legend
rlim([-30 5])

