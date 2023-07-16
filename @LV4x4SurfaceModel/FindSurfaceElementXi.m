function Xi = FindSurfaceElementXi(L,si,P)
% Calculate the surface element Xi(1), Xi(2) and Xi(3) for a given set of
% points in Cartesian coordinates.
%
%   Xi = L.FindSurfaceElementXi(si,P);
%
% Inputs:  - si is either 'endo' or 'epi'.
%          - P is Nx3 points in Cartesian coordinate
%
% Output:  - Xi is Nx4 matrix where each column is: 
%            [ Xi(1) Xi(2) Xi(3) element_ID ]
%            where Xi(1) is the element sampling in theta direction,
%                  Xi(2) is the element sampling in mu direction, and
%                  Xi(3) is the element sampling in lambda direction.
%
% The Xi output is similar with L.E matrix.
%
% Note: If points are not in the element ID, the element_ID value is -1.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

Xi = [];
if( isempty(L) ), return; end

% to switch endo/epi in the global parameter
if( strcmpi(si,'endo') )
    si = 2;
elseif( strcmpi(si,'epi') )
    si = 1;
else
    error('The argument si must be either ''endo'' or ''epi''.');
end

% get the portion of muParams & thetaParams for the given surface
i = (si-1)*20 + (1:20);
muParams = L.muParams(i);        % 20 mu parameters
thetaParams = L.thetaParams(i);  % 20 theta parameters

% each parameter contains 4 corners of the element, and there are 4
% elements for each dimension.
%
% The 4x4 element corners:
%  1     2     3     4     1
%  5     6     7     8     5
%  9    10    11    12     9
% 13    14    15    16    13
% 17    18    19    20    17
%

% so for each muParams or thetaParams vector, this the map
mapR =[  1  2  6  5  1; ... %  1
         2  3  7  6  2; ... %  2
         3  4  8  7  3; ... %  3
         4  1  5  8  4; ... %  4
         5  6 10  9  5; ... %  5
         6  7 11 10  6; ... %  6
         7  8 12 11  7; ... %  7
         8  5  9 12  8; ... %  8
         9 10 14 13  9; ... %  9
        10 11 15 14 10; ... % 10
        11 12 16 15 11; ... % 11
        12  9 13 16 12; ... % 12
        13 14 18 17 13; ... % 13
        14 15 19 18 14; ... % 14
        15 16 20 19 15; ... % 15
        16 13 17 20 16];... % 16

% apply mapping
muParams = muParams(mapR);
thetaParams = thetaParams(mapR);

% the last column is the wrapping
thetaParams([4 8 12 16],[2 3]) = 2*pi;

% convert P into prolate coordinate
PC = LV4x4SurfaceModel.Cart2Prolate(L.focalLength,P(:,1),P(:,2),P(:,3));

% we're gonna use the inpolygon method from Matlab to determine which
% element it is
N = size(P,1);
elmtIJ = zeros(16,N);
for i=1:16
    [IN ON] = inpolygon(PC(:,1),PC(:,2),thetaParams(i,:)',muParams(i,:)');
    elmtIJ(i,:) = IN | ON;
end

[~,elmtID] = max(elmtIJ,[],1);

% now find Xi on each element
A = thetaParams(elmtID,2) - thetaParams(elmtID,3);
a = muParams(elmtID,2) - muParams(elmtID,3);

B = thetaParams(elmtID,4) - thetaParams(elmtID,3);
b = muParams(elmtID,4) - muParams(elmtID,3);

C = thetaParams(elmtID,3) - thetaParams(elmtID,4) - thetaParams(elmtID,2) + thetaParams(elmtID,1);
c = muParams(elmtID,3) - muParams(elmtID,4) - muParams(elmtID,2) + muParams(elmtID,1);

D = PC(:,1) - thetaParams(elmtID,3);
d = PC(:,2) - muParams(elmtID,3);

alpha = A.*a - a.*C;
beta = d.*C - D.*c + A.*b - a.*B;
gamma = B.*d - b.*D;

% Solve for xi (mu)

elmtMu = -(gamma ./ beta);

iN0 = alpha~=0;
elmtMu(iN0) = -beta(iN0) + sqrt(beta(iN0).^2-4.*alpha(iN0).*gamma(iN0))/2 .* alpha(iN0);

% if mu outside [0-1]
iOut = iN0 & (elmtMu<0 | elmtMu>1);
elmtMu(iOut) = -beta(iOut) - sqrt(beta(iOut).^2-4.*alpha(iOut).*gamma(iOut))/2 .* alpha(iOut);

% Find xi (theta)
% Compensate for divide by zero
i0 = b+c.*elmtMu == 0;
elmtTheta(i0) = (D(i0)-A(i0).*elmtMu(i0)) ./ (B(i0)+C(i0).*elmtMu(i0));

in0 = b+c.*elmtMu ~= 0;
elmtTheta(in0)=(d(in0)-a(in0).*elmtMu(in0)) ./ (b(in0)+c(in0).*elmtMu(in0));

% define Xi
Xi(:,1) = elmtTheta;
Xi(:,2) = elmtMu;
% mind the switching
if( si==2 ), Xi(:,3) = 0; elseif( si==1 ), Xi(:,3)=1; else error('INTERNAL ERROR'); end
Xi(:,4) = elmtID;
