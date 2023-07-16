function BrowseSlices(L)
% Browse LV model throughout SA & LA slices interactively.
%
%   L.BrowseSlices;
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute 2011

if( isempty(L) ), return; end

% set figure
h.fig = figure;

%axprop = {'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','on',...
%    'XTick',[],'YTick',[],'ZTick',[],'Color','none'};

h.ax3d = axes('Position',[0.05 0.05 0.45 0.9]); 
axis(h.ax3d,'vis3d','off'); 
view(-79,-70);

% plot the LV and get the colors
h.hobj.L3d = plot(L,'with_points',false,'with_wireframe',false); hold on;
for i=1:numel(h.hobj.L3d.surf)
    h.surfcolor(i,:) = get(h.hobj.L3d.surf(i),'FaceColor');
end

% LV axes
h.hobj.ax3d(1) = plot3(xlim,[0 0],[0 0],'-k','LineWidth',1.5);
h.hobj.ax3d(2) = plot3([0 0],ylim,[0 0],'-k','LineWidth',1.5);
h.hobj.ax3d(3) = plot3([0 0],[0 0],zlim,'-k','LineWidth',1.5);

h.axla = axes('Position',[0.55 0.05 0.4 0.425],'Color','k');
h.laangle = -1;
h.hobj.lapt2 = NaN*ones(size(h.surfcolor));
h.hobj.lapt3 = h.hobj.lapt2;
hold(h.axla,'on');

h.axsa = axes('Position',[0.55 0.525 0.4 0.425],'Color','k');
h.sapos = -1;
h.hobj.sapt2 = NaN*ones(size(h.surfcolor));
h.hobj.sapt3 = h.hobj.sapt2;
hold(h.axsa,'on');

% set other props
h.L = L;
guidata(h.fig,h);

% set key interaction
set(h.fig,'KeyPressFcn',{@a_key_is_pressed_callback});

modify_slices(h,0,0);

% last touches
axis(h.axsa,'image');
xlim(h.axsa,get(h.axsa,'XLim') + [-5 5]);
ylim(h.axsa,get(h.axsa,'YLim') + [-5 5]);
axis(h.axla,'image','ij');
xlim(h.axla,get(h.axla,'XLim') + [-5 5]);
ylim(h.axla,get(h.axla,'YLim') + [-5 5]);

% set handle
h = guidata(h.fig);
h.hobj.saguide = plot(get(h.axla,'XLim'),[h.sapos h.sapos],'LineSTyle','-','Color','y','Parent',h.axla);
P = [cos(h.laangle) -sin(h.laangle); sin(h.laangle) cos(h.laangle)] * [2*get(h.axsa,'XLim'); 0 0];
h.hobj.laguide = plot(P(1,:),P(2,:),'LineSTyle','-','Color','y','Parent',h.axsa);
guidata(h.fig,h);

end

% --- MODIFY_SLICES
function modify_slices(h,sapos,laangle)

if( h.sapos ~= sapos )
    
    h.sapos = sapos;
    P = h.L.GetIntersectionWithPlane([h.sapos 0 0],[1 0 0]);
    
    for i=1:numel(P)
        if( ~ishandle(h.hobj.sapt2(i)) )
            if( ~isempty(P{i}) )
                h.hobj.sapt2(i) = plot([P{i}(:,2); P{i}(1,2)],[P{i}(:,3); P{i}(1,3)],...
                    'Color',h.surfcolor(i,:),'Marker','.','LineStyle','-','Parent',h.axsa); 
                h.hobj.sapt3(i) = plot3([P{i}(:,1); P{i}(1,1)],[P{i}(:,2); P{i}(1,2)],...
                    [P{i}(:,3); P{i}(1,3)],'LineStyle','-','Color',h.surfcolor(i,:),'Parent',h.ax3d,'LineWidth',1.5);
            end
        else
            if( isempty(P{i}) ), delete([h.hobj.sapt2(i) h.hobj.sapt3(i)]);
            else
                set(h.hobj.sapt2(i),'XData',[P{i}(:,2); P{i}(1,2)],'YData',[P{i}(:,3); P{i}(1,3)]);
                set(h.hobj.sapt3(i),'XData',[P{i}(:,1); P{i}(1,1)],'YData',[P{i}(:,2); P{i}(1,2)],'ZData',[P{i}(:,3); P{i}(1,3)]);
                
            end
            set(h.hobj.saguide,'YData',[h.sapos h.sapos],'XData',get(h.axla,'XLim'));
        end
    end
end

if( h.laangle ~= laangle )
    
    h.laangle = laangle;
    P = h.L.GetIntersectionWithPlane([0 0 0],[0 cos(h.laangle) sin(h.laangle)]);
    R = inv([1 0 0; 0 cos(h.laangle) sin(h.laangle); 0 -sin(h.laangle) cos(h.laangle)]'); % rotation along X axes
    
    for i=1:numel(P)
        Pr = P{i} * R';
        if( ~ishandle(h.hobj.lapt2(i)) )
            if( ~isempty(P{i}) )
                h.hobj.lapt2(i) = plot(Pr(:,3),Pr(:,1),...
                    'Color',h.surfcolor(i,:),'Marker','.','LineStyle','-','Parent',h.axla); 
                h.hobj.lapt3(i) = plot3(P{i}(:,1),P{i}(:,2),P{i}(:,3),...
                    'LineStyle','-','Color',h.surfcolor(i,:),'Parent',h.ax3d,'LineWidth',1.5);
            end
        else
            if( isempty(P{i}) ), delete([h.hobj.lapt2(i) h.hobj.lapt3(i)]);
            else
                set(h.hobj.lapt2(i),'XData',Pr(:,3),'YData',Pr(:,1));
                set(h.hobj.lapt3(i),'XData',P{i}(:,1),'YData',P{i}(:,2),'ZData',P{i}(:,3));
                
            end
            G = [cos(h.laangle) -sin(h.laangle); sin(h.laangle) cos(h.laangle)] * [2*get(h.axsa,'XLim'); 0 0];
            set(h.hobj.laguide,'XData',G(1,:),'YData',G(2,:));
        end
    end
end

% update guidata
guidata(h.fig,h);

end


% --- A_KEY_IS_PRESSED_callback
function a_key_is_pressed_callback(hobj,K)

h = guidata(hobj);

if( numel(K.Modifier)==1 && strcmpi(K.Modifier{:},'shift') )
    step = 5; 
else
    step = 1;
end

% check the key
if( strcmpi(K.Key,'downarrow') )
   
    
    % to the apex (positive x)
    new_sa = h.sapos + step;
    if( new_sa < max(get(h.ax3d,'XLim')) ), modify_slices(h,new_sa,h.laangle); end
    
elseif( strcmpi(K.Key,'uparrow') )
    
    % to the base (negative x)
    new_sa = h.sapos - step;
    if( new_sa > min(get(h.ax3d,'XLim')) ), modify_slices(h,new_sa,h.laangle); end
    
elseif( strcmpi(K.Key,'leftarrow') )
    
    new_la = h.laangle - step*(pi/90);
    modify_slices(h,h.sapos,new_la);
    
elseif( strcmpi(K.Key,'rightarrow') )

    new_la = h.laangle + step*(pi/90);
    modify_slices(h,h.sapos,new_la);
    
elseif( strcmpi(K.Key,'0') )
    
    modify_slices(h,0,0);
    
end

end