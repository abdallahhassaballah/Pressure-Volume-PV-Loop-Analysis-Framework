function TML = Cart2Prolate( focalLength, X, Y, Z )
% Convert rectangular Cartesian coordinates into prolate spheroidal
% coordinates.
%
%   TML = Cart2Prolate(focalLength,X,Y,Z);
%
% Output: TML(:,1) = lambda (radial) parameters
%         TML(:,2) = mu (azimuth) parameters
%         TML(:,3) = theta (circumferential) parameters
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

r1=sqrt( Y.^2 + Z.^2 + (X+focalLength).^2 );
r2=sqrt( Y.^2 + Z.^2 + (X-focalLength).^2 );

P.lambda=real(acosh((r1+r2)./(2*focalLength)));
P.mu=real(acos((r1-r2)/(2*focalLength)));
P.theta = atan2(Z,Y);

% compensate -pi..0 for atan2 (the bottom half)
idx = P.theta<0;
P.theta(idx) = P.theta(idx) + 2*pi;

% WARNING:
% Theta at apex (Z,Y) = (0,0) will not recover their true angle.

TML = [P.theta(:) P.mu(:) P.lambda(:)];