function P = ComputeBasisFunctions(elmts,nBasis,basisType)
% Calculate the basis functions from the LV model defined by 4x4
% surface sampling (theta & mu) and the 2 surfaces: endo and epi.
%
%   P = ComputeBasisFunctions(elmts,nBasis,basisType)
%
% Input:  - elmts is a Nx4 matrix of N sampling points, where columns are:
%             [ dim1_samples, dim2_samples, dim3_samples, element_number]
%         - nBasis is the number of basis
%         - basisType value is the following:
%             'trilinear'
%             'bicubic-linear'
% Output: - P is a sparse matrix that contains the basis functions for each
%           element node
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)
% Note:
%   This script is a vectorized version of calcP_LV4x4PS.m created by
%   Alistair Young and David Cumin (2005).

% get the elements
xi = elmts(:,1:3);
nSamples = size(elmts,1);

% process for each basis type
if( strcmpi('trilinear',basisType) )

    % trilinear
    B(:,:,1) = 1-xi;
    B(:,:,2) = xi;

    % create a grid of i=1:2, j=1:2 and k=1:2, and multiply B on these grids
    [i,j,k] = ndgrid(1:2,1:2,1:2);
    bvalues = squeeze(B(:,1,i(:)) .* B(:,2,j(:)) .* B(:,3,k(:)));

    % first sparse column is [1 2 ... ngrids] repeated nsamples times
    sparse_i = (1:nSamples)' * ones(1,numel(i));
    
    % second sparse column contains the element nodes (see calcP_LV4x4PS.m)
    sparse_j = (elmts(:,4)-1) * nBasis * ones(1,numel(i)) + ...
        ones(nSamples,1) * ((j(:)-1)*2 + i(:) + 4*(k(:)-1))';

elseif( strcmpi('bicubic-linear',basisType) )
    
    % bicubic-linear
    xip2 = xi.*xi;
    xip3 = xip2.*xi;
    B(:,:,1) = 1-3*xi+3*xip2-  xip3;
    B(:,:,2) =   3*xi-6*xip2+3*xip3;
    B(:,:,3) =        3*xip2-3*xip3;
    B(:,:,4) =                 xip3;
    B(:,3,1) = 1-xi(:,3);
    B(:,3,2) = xi(:,3);
    
    % create a grid of i=1:4, j=1:4 and k=1:2, and multiply B on these grids
    [i,j,k] = ndgrid(1:4,1:4,1:2);
    bvalues = squeeze(B(:,1,i(:)) .* B(:,2,j(:)) .* B(:,3,k(:)));

    % first sparse column is [1 2 ... ngrids] repeated nsamples times
    sparse_i = (1:nSamples)' * ones(1,numel(i));

    % second sparse column contains the element nodes (see calcP_LV4x4PS.m)
    sparse_j = (elmts(:,4)-1) * nBasis * ones(1,numel(i)) + ...
        ones(nSamples,1) * ((j(:)-1)*4 + i(:) + 16*(k(:)-1))';
    
else
    error('Unhandled basis type %s',basisType);
end

% construct the sparse matrix
P = sparse(sparse_i(:),sparse_j(:),bvalues(:));


