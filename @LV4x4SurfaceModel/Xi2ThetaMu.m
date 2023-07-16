function TM = Xi2ThetaMu( L, si, Xi, Ei )
% Calculate the corresponding (theta,mu) coordinate for the given Xi
% surface elements.
%
%   TM = L.Xi2ThetaMu(si,Xi,Ei);
%
% Inputs: - si is either using 'endo' or 'epi' elements definition.
%         - Xi is Nx2 element sampling. Columns are [theta_sampling mu_sampling].
%         - Ei is N-length of element IDs, from 1 to 16.
%
% Output: - TM is Nx2 corresponding [theta mu] coordinates of Xi.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

TM = [];
if( isempty(L) ), return; end

% check Ei
if( any(Ei(:)<1) || any(Ei(:))>16 )
    error('Invalid element ID. Elements must be between 1 to 16.');
end

% get the correct surface prolate spheroidal points
if( strcmpi(si,'endo') )
    tml = L.surfacePSPoints(1:(size(L.surfacePSPoints,1)/2),:);
elseif( strcmpi(si,'epi') )
    tml = L.surfacePSPoints((size(L.surfacePSPoints,1)/2+1):end,:);
else
    error('The argument si must be either ''endo'' or ''epi''.');
end

% we need to warp Xi one element by one element
uEi = unique(Ei);
TM = zeros(size(Xi));
for i=1:numel(uEi)
    
    idx = Ei==uEi(i);
    TM(idx,:) = ThinPlateSpline([0 0; 0 1; 1 0; 1 1],tml(L.nodes(:,uEi(i)),1:2),Xi(idx,:));
    
end
