%Filename: - FitGeomModel.m
% Author: David Cumin
% Created: 1 August 2005
% Last modified: 1 August 2005

%  contains: function that makes sure there are no large 'jumps'
%   in the data from 3/2 pi to 0 and converts 0 to pi.
%            
function out = PhaseCorrect(in)
% Correct large 'jumps' in the phase domain, where 2pi wrap to 0.
%
%   out = PhaseCorrect(in)
%
% Author: David Cumin (2005)
% Modified by Avan for better reading (2011)

% append 0
out = [in(:); 0];

% for all points
for i=1:size(in,1)
    
    % alter if the distance between the two successive points > pi
    if( out(i+1) - out(i) > pi )
        
        % add 2pi
        out(i) = out(i) + 2*pi;
        
    end
    
end

% remove 0
out = out(1:size(in,1));