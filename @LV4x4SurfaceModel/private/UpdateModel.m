function L = UpdateModel(L)
% Update the LV model. Note this is a private method for LV4x4SurfaceModel
% class.
%
%   L = UpdateModel(L)
%
% Author: Avan Suinesiputra - Auckland Bioengineering Institute (2011)

% calculate the element points (mind phase correction)
lambdas = L.basis.Lambda * L.G2E.Lambda * L.lambdaParams;
mus = L.basis.Mu * L.G2E.Mu * L.muParams;
thetas = L.basis.Theta * PhaseCorrect(L.G2E.Theta * L.thetaParams);

L.surfacePSPoints = [thetas mus lambdas];

% now prolate to euclidean conversion
L.surfacePoints = LV4x4SurfaceModel.Prolate2Cart(L.focalLength,lambdas,mus,thetas);
