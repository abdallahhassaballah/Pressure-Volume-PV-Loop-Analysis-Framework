function E = RectangularSampling(nElements, nSamples)
% Generate samples between elements on a 2D surface without duplicate
% samples. The sample values are between 0 and 1 for each element.
%
%   E = RectangularSampling(nElements,nSamples)
%
% Input:  - elements is 2x1 vector to define the number of elements per
%           dimension.
%         - nSamples is 2x1 vector of the number of sample values per dimension.
% Output: - E will be a matrix of N x 3, where N is the total number of
%           samples and each column is:
%           - col 1 = sample values for dimension 1
%           - col 2 = sample values for dimension 2
%           - col 3 = element numbers
%
% Example:
%   E = RectangularSampling([4 4],[3 3]) creates 4x4 finite element
%   model of a 2D surface with each element contains 3x3 samples.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

if( numel(unique(nSamples))>1 )
    warning('BUG found when number of samples are not equal!!!');
end

nTotalElmts = prod(nElements);
nTotalSamples = prod(nSamples);

%% GenerateSamplingPoints
numElemts = nTotalElmts * nTotalSamples;
[e1 e2] = meshgrid(linspace(0,1,nSamples(1)),linspace(0,1,nSamples(2)));

% columns of E = [e1, e2, eid]
E = zeros(numElemts,3);
E(:,1) = repmat(e1(:),size(E,1)/numel(e1),1);
E(:,2) = repmat(e2(:),size(E,1)/numel(e2),1);
E(:,3) = reshape(ones(numel(e1),1) * (1:nTotalElmts),[],1);

%% remove duplicate nodes
idx = 1:nTotalSamples;             % full sampling on (1,1)
for i=2:nElements(2)   % full sampling only for the 2st dimension
    j = (i-1)*nTotalSamples;
    idx = [idx j+1:j+nTotalSamples-nSamples(1)];
end
for i=nElements(2)+1:nTotalElmts
    j = (i-1)*nTotalSamples;
    if( mod((i-1),nElements(2))==0 ) % full sampling only for the 1st dimensioin
        idx = [idx j+setdiff(1:nTotalSamples,nSamples(2):nSamples(2):nTotalSamples)];
    else % no full sampling
        idx = [idx j+setdiff(1:nTotalSamples-nSamples(1),nSamples(2):nSamples(2):nTotalSamples)];
    end
end

% return
E = E(idx,:);
