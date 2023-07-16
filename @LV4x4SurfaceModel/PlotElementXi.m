function h = PlotElementXi(L,si,P,varargin)
% Plot each element Xi in the element surface space (theta-mu)
%
%   h = PlotElementXi(si,P,'opt1',val1,...);
%
% Input: - si is either 'endo' or 'epi'.
%        - P is Cartesian coordinates
% Output: h is a graphics object handle.
%
% Available options:
%   - Other Matlab's plot options.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% get the xi
Xi = L.FindSurfaceElementXi(si,P);

% to switch endo/epi in the global parameter
% 1. muParams & thetaParams are [EPI; ENDO]
% 2. sufacePSPoints & surfacePoitns are [ENDO; EPI]
if( strcmpi(si,'endo') )
    tml = L.surfacePSPoints(1:(size(L.surfacePSPoints,1)/2),:);
elseif( strcmpi(si,'epi') )
    tml = L.surfacePSPoints((size(L.surfacePSPoints,1)/2+1):end,:);
else
    error('The argument si must be either ''endo'' or ''epi''.');
end

% we need to wrap Xi one element by one element
Xw = [];
for i=1:16
    
    Ps = Xi(Xi(:,4)==i,1:2);
    Pt = thin_plate_spline([0 0 0 1 1 0 1 1]',reshape(tml(L.nodes(:,i),1:2)',[],1),Ps);
    
    Xw = [Xw; Pt];
    
end

% set hold
on_hold = ishold;
if( ~on_hold ), hold on; end

h = plot(Xw(:,1),Xw(:,2),varargin{:});

% revert back hold status
if( ~on_hold ), hold off; end