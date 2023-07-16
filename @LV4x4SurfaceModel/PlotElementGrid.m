function h = PlotElementGrid(L,si,varargin)
% Plot the element grid, where x-axis is theta elements and y-axis is the
% mu elements.
%
%   h = L.PlotElementGrid(si,'opt1',val1,'opt2',val2,...);
%
% Input: si is either 'endo' or 'epi'.
%
% Optional options are all options for the Matlab's plot command.
% Output is the graphic's handle object of the grid.
%
% IMPORTANT NOTE: 
% ---------------
% The element grid is defined by L.muParams and L.thetaParams.
% Therefore, it depends on these global parameters. The layout of these
% vectos are different than surface elements & points. The first part of
% these parameters are for EPI and then the second part for ENDO. Whilst
% the element samplings are first for ENDO and then EPI.
%
% Based from CIM model, the following settings are fixed: L.muParams are 40x1
% and L.thetaParams are 40x1 vectors. The first 20 is for EPI, the second
% 20 is for ENDO.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

h = [];
if( isempty(L) ), return; end

% to switch endo/epi in the global parameter
if( strcmpi(si,'endo') )
    si = 2;
elseif( strcmpi(si,'epi') )
    si = 1;
else
    error('The argument si must be either ''endo'' or ''epi''.');
end

% check L.muParams & L.thetaParams
if( ~isequal(size(L.muParams),[40 1]) || ~isequal(size(L.thetaParams),[40 1]) )
    error('Either muParams or thetaParams are not 40x1 vectors. See help.');
end
    

% each parameter contains 4 corners of the element, and there are 4
% elements for each dimension.
%
% The 4x4 element corners:
%  1     2     3     4     1
%  5     6     7     8     5
%  9    10    11    12     9
% 13    14    15    16    13
% 17    18    19    20    17
%
% so for each muParams or thetaParams vector, this the map
mapR =[  1  2  6  5  1; ... %  1
         2  3  7  6  2; ... %  2
         3  4  8  7  3; ... %  3
         4  1  5  8  4; ... %  4
         5  6 10  9  5; ... %  5
         6  7 11 10  6; ... %  6
         7  8 12 11  7; ... %  7
         8  5  9 12  8; ... %  8
         9 10 14 13  9; ... %  9
        10 11 15 14 10; ... % 10
        11 12 16 15 11; ... % 11
        12  9 13 16 12; ... % 12
        13 14 18 17 13; ... % 13
        14 15 19 18 14; ... % 14
        15 16 20 19 15; ... % 15
        16 13 17 20 16];... % 16
        
% get the portion of muParams & thetaParams for the given surface
i = (si-1)*20 + (1:20);
muParams = L.muParams(i);        % 20 mu parameters
thetaParams = L.thetaParams(i);  % 20 theta parameters

% apply mapping
muParams = muParams(mapR);
thetaParams = thetaParams(mapR);

% the last column is the wrapping
thetaParams([4 8 12 16],[2 3]) = 2*pi;        

% plot
h = plot(thetaParams',muParams',varargin{:});

set(gca,'XTick',[0 0.5*pi pi 1.5*pi 2*pi],'XTickLabel',{'0','90','180','270','360'});
xlabel('theta');
ylabel('mu');
xlim([0 2*pi]);