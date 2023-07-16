function V = ComputeRegionalVolumes(L)
% LVMM/ComputeRegionalVolumes: compute endo and epicardial regional AHA volumes for each frames.
%
%   V = L.ComputeRegionalVolumes;
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

V.endo = zeros(16,L.nframes);
V.epi = zeros(17,L.nframes);

if( L.nframes<1 )
    return;
end

% ED frame to get the G indices
Li = L.model(L.ed);
G = Li.GetAHAPoints;

V.endo = [];
V.epi = [];

% calculate volumes 1..nFrames
for fi=1:L.nframes
    Li.SetParams(L.focalLengths(fi),L.lambdas(:,fi),L.mus(:,fi),L.thetas(:,fi));
    Vi = Li.CalculateRegionalAHAVolumes(G);
    
    V.endo = [V.endo Vi.endo(:)];
    V.epi = [V.epi Vi.epi(:)];
end 