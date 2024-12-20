function [A,B]=bem2d(xyb,topology,k,betaP,varargin);

% [A,B]=bem2d(xyb,topology,k,betaP{,chiefpoints});
%
% Calculates the coefficient matrix for the 2D BEM formulation.
% It admits a plane with finite or infinite impedance.
%
% Input variables:
%   -xyb:     node positions, first column is the x-coordinate, second
%             column is y-coordinate, and third column is the body
%             number to which the node belongs to.
%   -topology:each row contains the node numbers (row number in xyb) of the
%             nodes in one element. The last column is the body number the
%             element belongs to.
%   -k:       wavenumber.
%   -betaP:   normalised admittance of the plane, at k. If its value is 0,
%             a rigid plane is considered (infinite impedance); if its value
%             is NaN (not-a-number), free-field is assumed.
%   -chiefpoints: Like 'xyb', but contains CHIEF points instead
%             One row for each chief point
%
% Output variable:
%   -A:       coefficient matrix for the pressure.
%   -B:       coefficient matrix for the normal velocity.
%
% It is possible to calculate for the interior or exterior domain, or a
% combination ob both. The sign of the body numbers indicates whether the
% interior or exterior domain to that body must be considered.

% Susana Quiros y Alpera, Vicente Cutanda 5-2001.
% Vicente Cutanda Henriquez 03-2011, version including CHIEF points

[M,ncolxyb] = size(xyb);
[N,nknel] = size(topology);
nknel=nknel-1;
if ncolxyb<3
   error('Error: The input arrays must include body numbers - Calculation aborted');
end

if nargin<5
   chiefpoints=[];
else
   chiefpoints=varargin{1};
end
nchiefp=size(chiefpoints,1);

NumBodies = max(abs(xyb(:,3)));
if NumBodies ~= max(abs(topology(:,nknel+1)));
   error('Error: The input arrays must include consistent body numbers - Calculation aborted');
end

% Check geometry files to find dummy elements and nodes
if isnan(betaP)
   operate='n';
else
   operate='y';
end
[xybnd,topologynd,xydum,topodum,xynodum]=nodummy(xyb,topology,operate);
lndu=size(xybnd,1);
ledu=size(topologynd,1);

% coefficient matrix calculation
CConst=2*pi*(1+sign(xybnd(:,3)))/2; % 2*pi or 0 for exterior/interior domain.
%CConst=[CConst ; zeros(nchiefp,1)];
A=zeros(lndu+nchiefp,lndu);
B=zeros(lndu+nchiefp,lndu);
for jj=1:lndu+nchiefp
    if jj <= lndu
        pxyb=xybnd(jj,:);
    else
        pxyb=chiefpoints(jj-lndu,:);
    end
    for nde=1:N
        elknxyb=xyb(topology(nde,1:nknel),:);
        [F1A,F1B,CK]=intF1(pxyb,elknxyb,k,betaP);
        if jj <= lndu
            CConst(jj)=CConst(jj)-CK*(1+xydum(jj));% Nodes in contact with the plane have different C.
        end
        if topodum(nde)~=0
            A(jj,topologynd(topodum(nde),1:nknel))=A(jj,topologynd(topodum(nde),1:nknel))+F1A(1:nknel);
            B(jj,topologynd(topodum(nde),1:nknel))=B(jj,topologynd(topodum(nde),1:nknel))+F1B(1:nknel);
        end
    end
    if jj <= lndu
        disp([' 2D BEM calculation, row ' num2str(jj) ' of ' num2str(lndu+nchiefp) ...
            ', C constant = ' num2str(CConst(jj)/pi) ' * pi'])
        A(jj,jj)=A(jj,jj)-CConst(jj);
    end
end
