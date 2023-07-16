function EstimateEDES(L)
% Estimate ED & ES frame by calculating frame index of the maximum and minimum EDV.
% This function will change the L.ed and L.es values.
%
%   L.EstimateEDES;
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( L.nframes==0 ), return; end

% compute volumes
edv = zeros(1,L.nframes);
for fi=1:L.nframes
    edv(fi) = L.model(fi).GetEndoVolume;
end

% update
[~,L.ed] = max(edv);
[~,L.es] = min(edv);