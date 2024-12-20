%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021

function [G,cutoff_f] = greens_function_duct(f, x, y, z, m_max, n_max, x0, y0, a, b)
    % ''' calculate Green's function for an semi-infinite rectangular duct
    % assuming coordinate origin at (0,0)
    % either frequency || (x,y,z) can be arrays
    % inputs:
    %     f     : array-like (or double), frequency in [Hz]
    %     x     : double (or array-like), x coordinates of receiver, in [m]
    %     y     : double (or array-like), y coordinates of receiver, in [m]
    %     z     : double (or array-like), z coordinates of receiver, in [m]
    %     x0    : double, x coordinate of source position, in [m]
    %     y0    : double, y coordinate of source position, in [m]
    %     a     : double, size of the duct along dimension x, in [m]
    %     b     : double, size of the duct along dimension y, in [m]
    %     m_max : maximum mode order along x dimension, non-negative integer
    %     n_max : maximum mode order along y dimension, non-negative integer
    % '''

    if (a<0) || (b<0)
        error('ValueError: duct dimensions be non-negative')
    elseif any(x<0) || any(x>a)
        error('ValueError: x needs to be between 0 <= x <= a')
    elseif any(y<0) || any(y>b)
        error('ValueError: y needs to be between 0 <= y <= b')
    elseif (x0<0) || (x0>a)
        error('ValueError: source position x0 needs to be between 0 <= x <= a')
    elseif (y0<0) || (y0>b)
        error('ValueError: source position y0 needs to be between 0 <= y <= b')
    elseif any(z<0)
        error('ValueError: z needs to be between 0 <= y')
    elseif (m_max<0) || (n_max<0)
        error('ValueError: orders m and n must be non-negative integers')
    end

    % calculate wavenumber 
    c = 343;
    k = 2*pi*f / c;
    
    G = 0; % initialize array
    S = a*b;
    cutoff_f=zeros(m_max+1,n_max+1);
    
    % calculate and sum contributions up to the (m_max, n_max)'th mode
    for mm = 0:m_max
        for nn = 0:n_max
            %% YOUR CODE GOES HERE
            % calculate cutoff frequency for this mode
             fcmn = cutoff_frequency(mm, nn, a, b);

            % calculate axial wavenumber for this mode
             k_zmn = axial_wavenumber(k,mm,nn,a,b);

            % calculate contribution of this mode
            % mode_contribution = ...

            %%
            
            Psi = rectangular_mode_shape(x, y, mm, nn, a, b);
            Psi0 = rectangular_mode_shape(x0, y0, mm, nn, a, b);
            mode_contribution = -1i/S * (Psi*Psi0)./k_zmn.*exp(-1i*k_zmn*z);
            G = G + mode_contribution; % accumulate contributions in G
            cutoff_f(mm+1,nn+1)=fcmn;
            fprintf('mode (%d,%d) out of (%d,%d), cutoff at %.1f Hz\n', ...
                mm, nn, m_max, n_max, fcmn);
        end
    end

    %% YOUR CODE GOES HERE
    %G = G*1i*omega*rho*Q; %scale result

    %%
end
