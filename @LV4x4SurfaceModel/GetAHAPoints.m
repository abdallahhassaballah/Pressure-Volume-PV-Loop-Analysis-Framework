function [R,mus,thetas] = GetAHAPoints(L,varargin)
% Get points for each AHA segments
%
%   R = L.GetAHAPoints;
%   [R,mus,thetas] = L.GetAHAPoints('opt1',val1,...);
%
% Outputs:
%   - R contains 17 cells of indices of each point belongs to which region.
%   - mus are struct of 'epi' and 'endo' that contains slice level
%     for apical tip (epi only), apex, mid and base regions.
%   - thetas contains angles for boundaries between regions 1-16.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

R = [];
mus = struct('epi',[],'endo',[]);
thetas = [];
if( isempty(L) ), return; end

% do not use transformation
Ttmp = L.T;
L.T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% options
opt.rv_inserts = [];

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option.'); 
    end
end

% get the levels first
[mus,lvl_idx] = L.GetAHASliceLevels;

% define the thetas
if( ~isempty(opt.rv_inserts) )
        
    xirv = L.FindSurfaceElementXi('epi',opt.rv_inserts);
    xirv_epi = L.Xi2ThetaMu('epi',xirv(:,1:2),xirv(:,4));
    xirv_endo = L.Xi2ThetaMu('endo',xirv(:,1:2),xirv(:,4));
    
    i_inferior = xirv(:,4) == 1 | xirv(:,4) == 5 | xirv(:,4) == 9 | xirv(:,4) == 13;
    i_anterior = xirv(:,4) == 4 | xirv(:,4) == 8 | xirv(:,4) == 12 | xirv(:,4) == 16;
    
    rv_inferior = mean([xirv_epi(i_inferior(:)); xirv_endo(i_inferior(:))]);
    rv_anterior = mean([xirv_epi(i_anterior(:)); xirv_endo(i_anterior(:))]);
    
    if( isnan(rv_inferior) ), rv_inferior = pi/3; end
    if( isnan(rv_anterior) ), rv_anterior = 5*pi/3; end

    thetas = [repmat([4*pi/3 rv_anterior 0 rv_inferior 2*pi/3 pi],1,2) 5*pi/4 7*pi/4 pi/4 3*pi/4];
    
else
    thetas = [repmat([4*pi/3 5*pi/3 0 pi/3 2*pi/3 pi],1,2) 5*pi/4 7*pi/4 pi/4 3*pi/4];
end

% select indices that belong to each regions
R = cell(1,17);

% the apical has been defined from lvl_idx
R{17} = lvl_idx.epi{1};

% NOTE: !!!!
% The LV model is flipped (180 deg rotation) along the Y-axes
% so: [1 2 3 4 5 6] -> [4 3 2 1 6 5]

% 1 - 6
th1_6 = [thetas([3 4 5 6 1 2]) 2*pi+0.1];
base_idx = [lvl_idx.epi{4}; lvl_idx.endo{3}];
for i=1:numel(th1_6)-1
    idx = find(L.surfacePSPoints(:,1)>=th1_6(i) & L.surfacePSPoints(:,1)<th1_6(i+1));
    R{i} = intersect(idx,base_idx);
end
R(1:6) = R([5 6 1 2 3 4]);

% 7 - 12
th7_12 = [thetas([9 10 11 12 7 8]) 2*pi+0.1];
mid_idx = [lvl_idx.epi{3}; lvl_idx.endo{2}];
for i=1:numel(th7_12)-1
    idx = find(L.surfacePSPoints(:,1)>=th7_12(i) & L.surfacePSPoints(:,1)<th7_12(i+1));
    R{i+6} = intersect(idx,mid_idx);
end
R(7:12) = R([11 12 7 8 9 10]);

% 13 - 16
th13_16 = [0 thetas([15 16 13 14]) 2*pi+0.1];
apex_idx = [lvl_idx.epi{2}; lvl_idx.endo{1}];
for i=1:numel(th13_16)-1
    idx = find(L.surfacePSPoints(:,1)>=th13_16(i) & L.surfacePSPoints(:,1)<th13_16(i+1));
    apex_tmp{i} = intersect(idx,apex_idx);
end
R(13:16) = [apex_tmp(4), {[apex_tmp{1}; apex_tmp{end}]}, apex_tmp(2), apex_tmp(3)];


% restore transform
L.T = Ttmp;

% make a grid of R
R = AHAGrids(L,R);