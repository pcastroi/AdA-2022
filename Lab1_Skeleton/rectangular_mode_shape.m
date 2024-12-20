%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021

function Psi = rectangular_mode_shape(x, y, m, n, a, b)
    % ''' calculate mode shape across rectangular surface
    % assuming origin at (0,0)
    % inputs
    %     x : array-like, x coordinates of N positions, in [m]
    %     y : array-like, y coordinates of N positions, in [m]
    %     m : mode order along x dimension, non-negative integer
    %     n : mode order along y dimension, non-negative integer
    %     a : float, size of the duct along dimension x, in [m]
    %     b : float, size of the duct along dimension y, in [m]
    % outputs
    %     Psi_{m,n}(x,y) : mode shape of (m,n)'th mode at positions (x,y)
    % '''

    % check inputs are valid
    if (a<0) || (b<0)
        error('ValueError: duct dimensions be non-negative')
    elseif any(x<0) || any(x>a)
        error('ValueError: x needs to be between 0 <= x <= a')
    elseif any(y<0) || any(y>b)
        error('ValueError: y needs to be between 0 <= y <= b')
    elseif any(m<0) || any(n<0) || any(mod(m,1) ~= 0) || any(mod(n,1) ~=0)
        error('ValueError: orders m and n must be non-negative integers')
    end

    %% YOUR CODE GOES HERE ...
     
    if m == 0
        epsilon_m = 1;
    else
        epsilon_m = 2;
    end
    
    if n == 0
        epsilon_n = 1;
    else
        epsilon_n = 2;
    end 
     
    Psi = sqrt(epsilon_m*epsilon_n)*cos(m*pi*x/a)*cos(n*pi*y/b);


    %%
end
