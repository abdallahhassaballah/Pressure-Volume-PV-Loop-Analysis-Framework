function V = ComputeVolumes(L)
% LVMM/ComputeVolumes: compute endo and epicardial volumes for each frames.
%
%   V = L.ComputeVolumes;
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

V.endo = zeros(1,L.nframes);
V.epi = zeros(1,L.nframes);

if( L.nframes<1 )
    return;
end

% first frame
Li = L.model(1);
V.endo(1) = Li.GetEndoVolume;
V.epi(1) = Li.GetEpiVolume;

% calculate volumes 2..nFrames

for fi=2:L.nframes
    Li.SetParams(L.focalLengths(fi),L.lambdas(:,fi),L.mus(:,fi),L.thetas(:,fi));
    V.endo(fi) = Li.GetEndoVolume;
    V.epi(fi) = Li.GetEpiVolume;
end

myo = V.epi - V.endo;
mass = mean(myo);
mass = mass *1.04;
