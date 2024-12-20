%% ADVANCED ACOUSTICS LAB EXERCISE 1 %%%%%%%%%%%%%%%%%%%%%%%
% -- Green's function in a semi-infinite rectangular duct --
% manha@dtu.dk, feb 10th, 2021

function k_zmn = axial_wavenumber(k,m,n,a,b)
    % ''' returns the axial wavenumber 'k_zmn' for the (m,n)'th mode excited with
    % the wavenumber k
    % inputs: 
    %     k : wavenumber of the excitation
    %     m : mode order along x dimension, non-negative integer
    %     n : mode order along y dimension, non-negative integer
    %     a : float, size of the duct along dimension x, in [m]
    %     b : float, size of the duct along dimension y, in [m]
    % outputs :
    %     k_zmn : axial wavenumber
    % '''

    if any(k<0)
        error('ValueError: wavenumber k must be non-negative')
    elseif (a<0) || (b<0)
        error('ValueError: duct dimensions be non-negative')
    elseif any(m<0) || any(n<0) || any(mod(m,1) ~= 0) || any(mod(n,1) ~=0)
        error('ValueError: orders m and n must be non-negative integers')
    end

    %% YOUR CODE GOES HERE ...
    k_zmn = -sqrt(k.^2-(m*pi/a)^2-(n*pi/b)^2);
%     if imag(k_zmn) ~= 0
%         k_zmn = -k_zmn;
%     end
    % ensure sign is correct for negative sqrt arguments
    

    %%
end
