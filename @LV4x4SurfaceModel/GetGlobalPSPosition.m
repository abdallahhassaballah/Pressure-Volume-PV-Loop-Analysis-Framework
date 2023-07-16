function TML = GetGlobalPSPosition(L)
% Calculate the three coordinates of (theta,mu,lambda) of each 134 global
% coordinate positions.
%
%   TML = L.GetGlobalPSPosition;
%
% The 134 positions follow the following model topology of the lambda's Bezier
% points:
%
% Epicardial surface indices of lambda:
%  1-------------2-----------3-----------4---------1
%  |  1  2  6  5 | x x 10  9 | x x 14 13 | x x x x |
%  |  3  4  8  7 | x x 12 11 | x x 16 15 | x x x x |
%  | 19 20 24 23 | x x 32 31 | x x 28 27 | x x x x |
%  | 17 18 22 21 | x x 26 25 | x x 30 29 | x x x x |
%  5-------------6-----------7-----------8---------5
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  | 35 36 40 39 | x x 44 43 | x x 48 47 | x x x x |
%  | 33 34 38 37 | x x 42 41 | x x 46 45 | x x x x |
%  9------------10----------11----------12---------9
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  | 51 52 56 55 | x x 60 59 | x x 64 63 | x x x x |
%  | 49 50 54 53 | x x 58 57 | x x 62 61 | x x x x |
% 13------------14----------15----------16--------13
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  |  x  x  x  x | x x  x  x | x x  x  x | x x x x |
%  | 66  x  x 67 | x x  x  x | x x  x  x | x x x x |
%  | 65  x  x  x | x x  x  x | x x  x  x | x x x x |
% 17------------18----------19----------20--------17
%
% The endocardial surface indices are epicardial points + 67.
% The total number of points in the output are 67*2 = 134; the same number
% of lambda parameters.
%
% Note that each point has 4 components: L, dL/dMu, dL/dT and d^2L/dMdT,
% where L, M, T are lambda, mu and theta.
% E.g., node 2 on the first element has L, dL/dMu, dL/dT and d^2L/dMdT at
% indices 5 7 6 and 8.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

%
TML = nan(134,3);
if( isempty(L) ), return; end

% get the nodal parameters
Nepi = L.GetNodalParameters('epi');
Nendo = L.GetNodalParameters('endo');

% store the node parameters
node_idx = [1 5 9 13 17 21 25 29 33 37 41 45 49 53 57 61 65];

TML(:,3) = L.lambdaParams;
TML(node_idx,1:2) = Nepi(1:17,1:2);
TML(67+node_idx,1:2) = Nendo(1:17,1:2);

%
% dL/dMu

dMuEpi = (Nepi(1:16,2) - Nepi(5:20,2))/3;
dMuEndo =(Nendo(1:16,2) - Nendo(5:20,2))/3;

pos_idx = [3 7 11 15];
node_idx = 1:4;
TML(pos_idx,1:2) = [Nepi(node_idx,1) Nepi(node_idx,2) - dMuEpi(node_idx)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1) Nendo(node_idx,2) - dMuEndo(node_idx)];

pos_idx = [19 23 31 27 35 39 43 47 51 55 59 63 66 67];
node_idx = 5:18;
TML(pos_idx,1:2) = [Nepi(node_idx,1) Nepi(node_idx,2) + dMuEpi(node_idx-4)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1) Nendo(node_idx,2) + dMuEndo(node_idx-4)];


%
% dL/dT

dTheta = pi/6;   % 0.5pi / 3

pos_idx = [2 18 34 50];
node_idx = [1 5 9 13];
TML(pos_idx,1:2) = [Nepi(node_idx,1) + dTheta Nepi(node_idx,2)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1) + dTheta Nendo(node_idx,2)];

pos_idx = [6 10 14 22 26 30 38 42 46 54 58 62];
node_idx = [2 3 4 6 7 8 10 11 12 14 15 16];
TML(pos_idx,1:2) = [Nepi(node_idx,1) - dTheta Nepi(node_idx,2)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1) - dTheta Nendo(node_idx,2)];

%
% d^2L/dMudT

pos_idx = 4;
node_idx = 1;
TML(pos_idx,1:2) = [Nepi(node_idx,1)-dTheta Nepi(node_idx,2)-dMuEpi(node_idx)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1)+dTheta Nendo(node_idx,2)-dMuEndo(node_idx)];

pos_idx = [20 36 52];
node_idx = [5 9 13];
TML(pos_idx,1:2) = [Nepi(node_idx,1)+dTheta Nepi(node_idx,2)+dMuEpi(node_idx)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1)+dTheta Nendo(node_idx,2)+dMuEndo(node_idx)];

pos_idx = [8 12 16];
node_idx = [2 3 4];
TML(pos_idx,1:2) = [Nepi(node_idx,1)+dTheta Nepi(node_idx,2)-dMuEpi(node_idx)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1)+dTheta Nendo(node_idx,2)-dMuEndo(node_idx)];

pos_idx = [24 32 28 40 44 48 56 60 64];
node_idx = [6 7 8 10 11 12 14 15 16];
TML(pos_idx,1:2) = [Nepi(node_idx,1)-dTheta Nepi(node_idx,2)+dMuEpi(node_idx)];
TML(67+pos_idx,1:2) = [Nendo(node_idx,1)-dTheta Nendo(node_idx,2)+dMuEndo(node_idx)];
