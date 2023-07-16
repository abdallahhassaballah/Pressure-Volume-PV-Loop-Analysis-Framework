function idx = ShiftEDToFirst(L)
% Shift ED frame position to the first frame
%
%   shifted_idx = L.ShiftEDToFirst;
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

idx = 1:L.nframes;
if( L.nframes==0 ), return; end
if( L.ed==1 ), return; end

if( ~isempty(fieldnames(L.data) ) )
    fprintf(2,'Warning: shifting ED to first might invalidate external data.\n');
    ys = input('Are you sure to continue (y/n)? ','s');
    if( ~strcmpi(ys,'y') ), return; end
end

% re-indexing
idx = [L.ed:L.nframes 1:(L.ed-1)];

L.focalLengths = L.focalLengths(idx);
L.lambdas = L.lambdas(:,idx);
L.thetas = L.thetas(:,idx);
L.mus = L.mus(:,idx);
if( numel(L.phase_data)==L.nframes )
    L.phase_data = L.phase_data(idx);
end

% update ed & es
if( L.ed < L.es )
    L.es = L.es - (L.ed-1);
else
    L.es = L.es + (L.nframes-L.ed+1);
end

L.ed = 1;
