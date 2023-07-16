function L = UpdateSampling(L)
% Update element sampling. This is called usually because of changes in the
% number of samples / surfaces properties.
%
%   L = L.UpdateSampling;
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

% create 4x4 elements on a surface
E = LV4x4SurfaceModel.RectangularSampling([4 4],L.nSamples);

% create surface sample values
surfSamples = ones(size(E,1),1) * linspace(0,1,L.nSurfaces);

% repeat E nSurfaces times
L.E = repmat(E,L.nSurfaces,1);

% swap column 3 & 4
L.E = [L.E(:,1:2) surfSamples(:) L.E(:,3)];    

% update basis functions
L.basis.Lambda = L.ComputeBasisFunctions(L.E,32,'bicubic-linear');    % lambda
L.basis.Mu = L.ComputeBasisFunctions(L.E,8,'trilinear');          % mu
L.basis.Theta = L.ComputeBasisFunctions(L.E,8,'trilinear');          % theta

% update map
[L.map, L.nodes] = CreateSamplingMap(L);

% update model
if( ~isempty(L) ), L = L.UpdateModel; end