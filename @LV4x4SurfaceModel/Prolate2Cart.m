function XYZ = Prolate2Cart( focalLength, lambda, mu, theta )
% Convert prolate spheroidal coordinates into rectangular Cartesian
% coordinate.
%
%   XYZ = Prolate2Cart( focalLength, lambda, mu, theta )
%
% Output is Nx3 matrix of the cartesian coordinates, where N is the number
% of elements in lambda, theta and mu vectors. The focalLength must be a
% single number.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

XYZ(:,1) = focalLength .* cosh(lambda(:)) .* cos(mu(:));
XYZ(:,2) = focalLength .* sinh(lambda(:)) .* sin(mu(:)) .* cos(theta(:));
XYZ(:,3) = focalLength .* sinh(lambda(:)) .* sin(mu(:)) .* sin(theta(:));
