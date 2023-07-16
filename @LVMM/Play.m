function Play(L,varargin)
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

% check if it is empty
if( L.nframes==0 ), return; end

% default options
opts.axes = [];     % where to plot
opts.start = true;  % start immediately

% get options
for i=1:2:length(varargin)
    if( isfield(opts,varargin{i}) ), opts.(varargin{i}) = varargin{i+1};
    else
        error('Unknown option.');
    end
end

% set axes
 if( isempty(opts.axes) )
     opts.axes = axes('Visible','on');
end

% create LV models
for i=1:L.nframes
    h.models(i) = L.model(i);
end

% setup guidata
h.fig = gcf;
h.axes = opts.axes;
h.surf = [];
h.cp = L.ed;
h.onplay = true;
h.timer = timer('TimerFcn',@timer_callback,'Period',1/10,'BusyMode','queue',...
    'ExecutionMode','fixedDelay','UserData',h.fig);

% store guidata
guidata(h.fig,h);

disp('Stop and play with ''P'' key.');

% plot model
UpdateModel(h.fig);

% key press
set(h.fig,'KeypressFcn',@keypress_callback);

% start play directly
if( opts.start ) start(h.timer); end

end


% -- Update the model
function UpdateModel(hobj)

% get the models data
h = guidata(hobj);

% new plot or update
if( isempty(h.surf) )
    h.surf = h.models(h.cp).PlotSurface('Parent',h.axes);
    axis(h.axes,'image');
    view(h.axes,-700,-700);
    guidata(h.fig,h);
else
    % change
    P = h.models(h.cp).PlotSurface('noplot',true);

    for i=1:numel(P)
    set(h.surf(i),'Vertices',P(i).FV.Vertices,'Faces',P(i).FV.Faces,...
        'FaceVertexCData',P(i).FV.FaceVertexCData);
    end
end

end


% --- KEYPRESS
function keypress_callback(hobject,evtdata)

start_the_motion = false;

h = guidata(hobject);
if( lower(evtdata.Character) == 'p' )
    if( h.onplay )
        disp('stopping motion.');
    else
        start_the_motion = true;
        disp('playing motion.');
    end
    h.onplay = ~h.onplay;
end
guidata(h.fig,h);

if( start_the_motion ), start(h.timer); end

end



% --- TIMER CALLBACK 
function timer_callback(tmrObj,~)

h = guidata(get(tmrObj,'UserData'));
if( h.onplay )
    h.cp = mod(h.cp,numel(h.models))+1;
    guidata(h.fig,h);

    UpdateModel(h.fig);
else
    stop(tmrObj);
end

end
