function P = GetSurfacePointsFromElements( L, E)
% Get surface points given new element E
%
%   P = L.GetSurfacePointsFromElements(E);
%
% Input:  - E is Nx4 elements, the same format as L.E
% Output: - P is Nx3 Cartesian coordinates of the surface points of E.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( isempty(L) )
    P = [];
    return;
end

% create basis
BLambda = L.ComputeBasisFunctions(E,32,'bicubic-linear');    % lambda
BMu = L.ComputeBasisFunctions(E,8,'trilinear');          % mu
BTheta = L.ComputeBasisFunctions(E,8,'trilinear');          % theta

% create P
lambdas = BLambda * L.G2E.Lambda * L.lambdaParams;
mus = BMu * L.G2E.Mu * L.muParams;
thetas = BTheta * PhaseCorrect(L.G2E.Theta * L.thetaParams);

% now prolate to euclidean conversion
P = LV4x4SurfaceModel.Prolate2Cart(L.focalLength,lambdas,mus,thetas);

% apply transformation
P = [P ones(size(P,1),1)] * L.T;
