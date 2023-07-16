function view(L)
% View LV4x4SurfaceModel object.
% 
%   view(L);
%
% This function is similar to plot but with extra interactive and menu features.
% View will create new figure and you can only view one object at a time.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( isempty(L) ), return; end

% create new figure
h.fig = figure('Color','w','NumberTitle','off','MenuBar','none','Name','');
h.ax = axes('Visible','off','DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual',...
    'PlotBoxAspectRatio',[1.1353 1 1.1508],'PlotBoxAspectRatioMode','manual');

% plot
hp = plot(L);
h.surf = hp.surf;
h.wire = [hp.wire.circ; hp.wire.long];
h.nodes = hp.wire.nodes;
h.points = hp.points;
h.show_endo_epi = [1 1];

% set toolbar
ht = uitoolbar('Parent',h.fig);

% set view menu
hv = uimenu(h.fig,'Label','View');
uimenu(hv,'Label','Surface','Checked','on','Callback',{@show_callback,'surf'});
uimenu(hv,'Label','Wireframe','Checked','on','Callback',{@show_callback,'wire'});
uimenu(hv,'Label','Surface points','Checked','on','Callback',{@show_callback,'points'});
uimenu(hv,'Label','Nodal points','Checked','on','Callback',{@show_callback,'nodes'});
uimenu(hv,'Label','Epicardium','Checked','on','Separator','on','Callback',{@show_endoepi_callback,2});
uimenu(hv,'Label','Endocardium','Checked','on','Callback',{@show_endoepi_callback,1});

% store guidata
guidata(h.fig,h);

end


% --- SHOW CALLBACK
function show_callback(hobj,~,x)

ischk = strcmpi(get(hobj,'Checked'),'on');
if( ischk ), set(hobj,'Checked','off');
else set(hobj,'Checked','on'); end

% get data
h = guidata(hobj);

% set only if endo/epi is set active
hx = h.(x);
if( strcmpi(get(hobj,'Checked'),'off') )
    set(hx,'Visible','off');
else
    set(hx(:,h.show_endo_epi==1),'Visible','on');
end

end

% --- SHOW ENDO EPI CALLBACK
function show_endoepi_callback(hobj,~,ix)

% get data
h = guidata(hobj);

if( strcmpi(get(hobj,'Checked'),'on') )
    h.show_endo_epi(ix) = 0;
    newval = 'off';
else
    h.show_endo_epi(ix) = 1;
    newval = 'on';
end
set(hobj,'Checked',newval);

objs = {'surf','wire','nodes','points'};
for i=1:numel(objs)
    hx = h.(objs{i});
    if( strcmpi(get(hx(:,ix),'Visible'),'on') )
        set(hx(:,ix),'Visible',newval);
    end
end

% store
guidata(h.fig,h);

end