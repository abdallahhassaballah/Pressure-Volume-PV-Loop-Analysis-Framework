function PlayAHA(L,varargin)
% Play LV motion movie with AHA regions shown and fixed orientation (flipped XY)
%
%   L.PlayAHA;
%
% Start & stop movie: 'P' key
% Next frame (will stop the movie): right arrow
% Preview frame (will stop the movie): left arrow
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( isempty(L) ), return; end


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
    opts.axes = axes('Visible','off');
end

% store model
h.Li = L.model(1);
h.Lm = L;
h.Li.T = [0 -1 0; -1 0 0; 0 0 1; 0 0 0];

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

h.Li.SetParams(h.Lm.focalLengths(h.cp),h.Lm.lambdas(:,h.cp),h.Lm.mus(:,h.cp),h.Lm.thetas(:,h.cp));

% new plot or update
if( isempty(h.surf) )
    h.surf = h.Li.PlotSurface('Parent',h.axes);
    axis(h.axes,'image');
    view(h.axes,-70,-70);
    guidata(h.fig,h);
else
    % change
    P = h.Li.PlotSurface('noplot',true);

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
    h.cp = mod(h.cp,h.Lm.nframes)+1;
    guidata(h.fig,h);

    UpdateModel(h.fig);
else
    stop(tmrObj);
end

end
