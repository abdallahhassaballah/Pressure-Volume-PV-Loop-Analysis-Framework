function [N,idx] = GetNodalParameters(L,si)
% Get the global prolate spheroidal parameters on the all the 4x4 nodes.
%
%   [N,idx] = L.GetNodalParameters(si)
%
% Input: si is the surface, either 'epi' or 'endo'.
% Output:
%   - N is 20x3 matrix, where each columns are global theta, mu and lambda
%     parameters. Each rows follow the surface map (from base to apex):
%
%       1 ----  2 ----  3 ----  4 ----  1      <-- base
%       |       |       |       |       |
%       5 ----  6 ----  7 ----  8 ----  5
%       |       |       |       |       |
%       9 ---- 10 ---- 11 ---- 12 ----  9
%       |       |       |       |       |
%      13 ---- 14 ---- 15 ---- 16 ---- 13
%       |       |       |       |       |
%      17 ---- 18 ---- 19 ---- 20 ---- 17      <-- apex
%
%       0      p/2     pi     3*pi/2   2*pi
%
%   - idx is 20x3 indices of the global parameters
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

N = zeros(20,3);
if( isempty(L) ), return; end

% get the surface offset q
if( strcmpi(si,'epi') ), qL = 0; qMT = 0;
elseif( strcmpi(si,'endo') ), qL = 67; qMT = 20;
else error('Possible value for ''si'' are ''epi'' or ''endo''.');
end

% indices
mtidx = qMT + (1:20)';
lidx = qL + [1:4:61 repmat(65,1,4)]';

% get the value
N = [L.thetaParams(mtidx) L.muParams(mtidx) L.lambdaParams(lidx)];
idx = [mtidx mtidx lidx];
