function  h = hankelsph(m,z)
    % Hankel function of second kind
    H = besselh(m+1/2,2,z);
    
    % cylindrical Hankel function of second kind
    h = sqrt(pi./(2.*z)).*H;
end