%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021

function cutoff_frequency = cutoff_frequency(m, n, a, b)
    % ''' returns the cutoff frequency in Hz for the (m,n)'th mode
    %     m : mode order along x dimension, non-negative integer
    %     n : mode order along y dimension, non-negative integer
    %     a : float, size of the duct along dimension x, in [m]
    %     b : float, size of the duct along dimension y, in [m]
    % '''

    if (a<0) || (b<0)
        error('ValueError: duct dimensions be non-negative')
    elseif any(m<0) || any(n<0) || any(mod(m,1) ~= 0) || any(mod(n,1) ~=0)
        error('ValueError: orders m and n must be non-negative integers')
    end

    %% YOUR CODE GOES HERE ...
    c = 343;
    cutoff_frequency = c/2*sqrt((m/a)^2+(n/b)^2);
    %%
end
