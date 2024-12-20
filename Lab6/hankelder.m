function dh = hankelder(m,z)
    dh = 1./(2*m+1).*(m.*hankelsph(m-1,z)-(m+1).*hankelsph(m+1,z));
end