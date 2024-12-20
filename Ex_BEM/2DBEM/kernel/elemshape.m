function [psi, xq, yq, nx, ny]=elemshape(elknxyb,IP)

% [psi, xq, yq, nx, ny]=elemshape(elknxyzb,IP)
%
% Calculate shape functions, global coordinates and normal vector at
% a given set of local points 'IP', where 'e' is the local
% coordinate in an isoparametric element, either linear or quadratic.
% Version for 2D BEM.
%
%  Input:
%    IP:       real vector containing the local coordinate for the 
%             integration points.
%    elknxyb: real matrix, one row for each node in the element
%             each row contains (x,y,body) for the node.
%
%  Output:
%    psi:     Real matrix containing one shape function in each
%             column (no. of columns=no. of nodes in an element).
%    xq:      Real column vector containing the global x-coordinates
%             of the element for each value in IP.
%    yq:      Real column vector containing the global y-coordinates
%             of the element for each value in IP.
%    nx,ny:   The column vectors containing respectively the global
%             x and y-coordinates of the normal vector.
%             Note that the normal vector does not have unit length.
%
% References: -Peter M. Juhl: "The Boundary Element Method for Sound
%             Field Calculations", Report No. 55, DTU 1993. Section 4.7.
%             -O. C. Zienkiewicz, R. L. Taylor: "The Finite Element Method"
%             4th Ed., Volume 1, Section 8.5.

% Peter M. Juhl, Susana Quiros y Alpera, Vicente Cutanda 2001

nknel=size(elknxyb,1);

if nknel==2
   % Linear Shape functions:
   psi=[0.5*(1-IP)  0.5*(1+IP)];
   % Linear Shape function derivatives
   dNde=[-0.5*ones(size(IP))  0.5*ones(size(IP))];
elseif nknel==3
   % Quadratic Shape functions:
   psi=[0.5*IP.*(IP-1)  1-IP.^2  0.5*IP.*(IP+1)];
   % Quadratic Shape function derivatives
   dNde=[IP-0.5  -2*IP  IP+0.5];
else
   error('Only linear & quadratic shape functions are implemented');
end

% Global coordinates of the integration points.
xq=psi*elknxyb(:,1);
yq=psi*elknxyb(:,2);

% Elements of the jacobian.
dxde=dNde*elknxyb(:,1);
dyde=dNde*elknxyb(:,2);

% Normal vector at integration points.
IntExt=sign(elknxyb(1,3)); % Change sign if the domain is interior.
nx=-dyde*IntExt;
ny=dxde*IntExt;
