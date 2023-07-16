function MyPlay(L,varargin)
% Play the LV motion.
%
%   L.Play('opt1',val1,...);
%
% Available options:
%   - 'axes', axes_handle. Default is the current axes.
%   - 'start', true|false. Default is true.
%     Start motion immediately.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

figure;
for i=1:L.nframes
    h.models(i) = L.model(i);
    subplot(5,5,i), h.models(i).PlotSurface;
    title (i);
end
%L.model(1).PlotSurface

end
