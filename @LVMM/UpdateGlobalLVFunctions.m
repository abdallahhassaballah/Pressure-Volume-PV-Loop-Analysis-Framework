function UpdateGlobalLVFunctions(L)
% Compute global LV functions: EDV, ESV, MASS, SV, and EF.
%
%   L.UpdateGlobalLVFunctions;
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( L.ed>=1 && L.es>=1 && L.ed<=L.nframes && L.es<=L.nframes )
    
    Led = L.model(L.ed);
    Les = L.model(L.es);
    V = L.ComputeVolumes;
    myo = mean((V.epi - V.endo));
    L.EDV = Led.GetEndoVolume;
    L.ESV = Les.GetEndoVolume;
    L.SV = L.EDV-L.ESV;
    L.EF = L.SV ./ L.EDV * 100; 
   % L.MASS = 1.04 * (Led.GetEpiVolume - L.EDV);
    L.MASS= 1.04*myo;
else
    
    L.EDV = NaN;
    L.ESV = NaN;
    L.SV = NaN;
    L.EF = NaN;
    L.MASS = NaN;
    
end