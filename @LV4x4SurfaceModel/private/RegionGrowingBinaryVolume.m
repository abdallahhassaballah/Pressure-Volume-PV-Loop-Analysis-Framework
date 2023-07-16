function [B,Xi,Yi,Zi] = RegionGrowingBinaryVolume(varargin)
% Generate a binary volume based on the seed growing method.
%
%   [B,Xi,Yi,Zi] = RegionGrowingBinaryVolume('opt1',val1,'opt2',val2,...);
%
% Available options:
%   - 'xlim' | 'ylim' | 'zlim', [min max].
%     Define limits of each axis. Default is [-10 10] for all x, y, and z
%     axes.
%   - 'incr', [xincr yincr zincr].
%     Define point increments on each axis. Default is [1 1 1].
%   - 'seed', seed_point.
%     Define the starting point of the growing iteration. Default is at the
%     origin, or [0 0 0].
%   - 'eval_fun', a function handle.
%     Define a function to evaluate points whether they are included
%     (binary value is 1) or exclude (=0). Default is just an identity
%     function.
%     This function takes 3 arguments: Xi, Yi and Zi, which are the three
%     coordinate values. The output of this function is a logical vector
%     that indicates whether each row of the input arguments are either
%     inside (true) or outside (false).
%
%     E.g., to create points inside 5 units radius, define eval_fun as:
%       f = @(xi,yi,zi) sqrt(xi*xi+yi*yi+zi*zi)<=5;
%
%     See doc function_handle.
%
% The output B is a three-dimensional points define by the limits and
% the increments where each element contains 0 (excluded) or 1 (included).
%
% The mesh grids are defined in Xi, Yi and Zi.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

% default options
opt.xlim = [-10 10];
opt.ylim = [-10 10];
opt.zlim = [-10 10];
opt.incr = [1 1 1];
opt.seed = [0 0 0];
opt.eval_fun = @(xi,~,~) ones(size(xi))==1;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

% create grids
[Xi,Yi,Zi] = ndgrid(...
    opt.xlim(1):opt.incr(1):opt.xlim(2), ...
    opt.ylim(1):opt.incr(2):opt.ylim(2), ...
    opt.zlim(1):opt.incr(3):opt.zlim(2));

% create visited points
B = zeros(size(Xi));

% index function
fidx = @(i,j,k) i.*size(B,1).*size(B,2) + j.*size(B,1) + k;

% neighbouring grids
[ii,jj,kk] = ndgrid([-1 0 1]*opt.incr(1),[-1 0 1]*opt.incr(2),[-1 0 1]*opt.incr(3));
neighbors = fidx(ii(:),jj(:),kk(:));

% range & origin
range = [1 numel(B)];
orig = find(Xi==opt.seed(1) & Yi==opt.seed(2) & Zi==opt.seed(3));
if( isempty(orig) ), error('Origin must be included in the grid.'); end

% front points
fronts = orig + fidx(opt.seed(1),opt.seed(2),opt.seed(3));
    
% iterations
while ~isempty(fronts)
    
    % visited fronts
    B(fronts) = 1;

    % get new fronts
    newf = unique(fronts(:) * ones(1,numel(neighbors)) + ones(numel(fronts),1) *  neighbors','rows');
    newf = newf(newf>=range(1) & newf<=range(2));
    fronts = newf(~B(newf));
    
    % apply evaluate function
    iEval = opt.eval_fun(Xi(fronts),Yi(fronts),Zi(fronts));
    fronts = fronts(iEval);
     
    % message
    fprintf(1,'RegionGrowingBinaryVolume: Added %d points.\n',numel(fronts));
    
end 

