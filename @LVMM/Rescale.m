function Rescale(L,new_focalLengths)
% Isotropic scale
%
%   L.Rescale(new_focalLengths);
%
% If numel(new_focalLengths)==1, then all models will be rescaled with one focalLength.
%
% Author: Avan Suinesiaputra - University of Auckland 2012

if( numel(new_focalLengths)==1 ), new_focalLengths = new_focalLengths*ones(size(L.focalLengths)); end
if( ~isequal(size(new_focalLengths),size(L.focalLengths)) )
    error('The number of new focalLengths does not match.');
end

L.focalLengths = new_focalLengths;