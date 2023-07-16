function L = CreateEndoEpiSurfacesByTML( focalLength, TML, varargin )
% Create 2 surfaces LV object by giving TML (theta-mu-lambda) matrix.
%
%   L = LV4x4SurfaceModel.CreateEndoEpiSurfacesByTML( focalLength, TML );
%   L = LV4x4SurfaceModel.CreateEndoEpiSurfacesByTML( focalLength, TML, ... );
%
% Input:  - focalLength is a scalar to define the model's focal length.
%         - TML is 134x3 matrix where columns are theta, mu and lambda
%           global parameters.
%
% Optional arguments: see LV4x4SurfaceModel.Create function.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

lambdaParams = TML(:,3);
muParams = [TML([1:4:61 65 65 65 65],2); TML(67+[1:4:61 65 65 65 65],2)];
thetaParams = [TML([1:4:61 49:4:61],1); TML(67+[1:4:61 49:4:61],1)];

L = LV4x4SurfaceModel.Create(focalLength,lambdaParams,muParams,thetaParams,varargin{:});
