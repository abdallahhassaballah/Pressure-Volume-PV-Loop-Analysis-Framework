function [SW, PE, PVA, eta,EEV, CO, Ea,Ees,Coupling_ratio,ESP,MEP] = Extract_PV_Loop_Paramaters(P,V,ESV,SV,HR)

Conversion_Factor = 0.000133322368; % mmHg*mL ==> J

% 1. Mitral Valve Opening (MVO): bottom left 
% This represents the start of diastole, where the left ventricle is relaxed 
% and filling with blood from the left atrium. At this point, the pressure 
% is relatively low, but the volume is increasing.

% 2. Mitral Valve Closure (MVC): bottom right
% As the left ventricle finishes filling with blood, 
% the pressure starts to increase. When the pressure in the left ventricle 
% exceeds that in the left atrium, the mitral valve closes, 
% marking the end of diastolic filling. This is the onset of isovolumetric contraction.

% 3. Aortic Valve Opening (AVO): top right
% When the pressure in the left ventricle 
% exceeds the pressure in the aorta, the aortic valve opens, allowing blood 
% to be ejected from the left ventricle into the aorta. 
% This is the start of the systolic ejection phase.

% 4. Aortic Valve Closure (AVC): top left
% As the pressure in the left ventricle 
% starts to drop towards the end of systole, when it becomes less than the 
% pressure in the aorta, the aortic valve closes. 
% This marks the end of systole and the beginning of isovolumetric relaxation.
% Calculate the Stroke Work (SW) as the area enclosed by the PV loop


% Find the four corners
min_volume = min(V);
max_volume = max(V);
min_pressure = min(P);
max_pressure = max(P);

% Find the indices of the data points closest to the four corners
[~, MVO_index] = min(abs(V-min_volume) + abs(P-min_pressure));
[~, AVC_index] = min(abs(V-min_volume) + abs(P-max_pressure));
[~, MVC_index] = min(abs(V-max_volume) + abs(P-min_pressure));
[~, AVO_index] = min(abs(V-max_volume) + abs(P-max_pressure));


SW = abs(trapz(V, P))* Conversion_Factor;


% Initialize the coordinates of the three points
x1 = 0; y1 = 0; % P(0), V(o)
x2 = V(AVC_index); y2 = P(AVC_index); % V(AVC_index), P(AVC_index)
x3 = V(MVO_index); y3 = P(MVO_index); % V(MVO_index), P(MVO_index)

% Calculate the area of the triangle
area = 0.5 * abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));

triangle_X = [x1, x2, x3, x1];
triangle_Y = [y1, y2, y3, y1];


% Delete the overlapping area between PE_org and SW
[x_int, y_int] = polybool('subtraction', triangle_X, triangle_Y, V, P);
PE = polyarea(x_int, y_int)* Conversion_Factor;

ESP = P(AVC_index);


% Calculate PV area
PVA = SW + PE;

% Calculate ventricular efficiency
eta = SW / PVA;

%The energy per ejected volume (EEV)
EEV = PVA/SV;

% The mean external power (MEP) that the LV delivers
MEP = SW*(HR/60);


CO = SV*HR;
Ea = ESP/SV;
Ees = ESP/ESV;

Coupling_ratio = SV/ESV;

%{
figure; hold on
plot(V, P, 'k') % plot the PV loop
scatter(V(MVO_index),P(MVO_index), 'r', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5) % plot the bottom left corner
scatter(V(MVC_index),P(MVC_index), 'g', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5) % plot the bottom right corner
scatter(V(AVO_index),P(AVO_index), 'b', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5) % plot the top right corner
scatter(V(AVC_index),P(AVC_index), 'y', 'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5) % plot the top left corner
%g = fill(triangle_X, triangle_Y, 'k', 'FaceAlpha', 0.1);
g = fill(x_int, y_int, 'k', 'FaceAlpha', 0.3);
h = fill(V, P, 'r', 'FaceAlpha', 0.3); 
legend('PV loop', 'Mitral Valve Opening', 'Mitral Valve Closure', 'Aortic Valve Opening ', 'Aortic Valve Closure', 'PE','SW')
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
title('Pressure-Volume Relationships')


ESV = min(V);
EDV = max(V);
ESP = P(AVC_index);
EDP = P(MVC_index);

%}

end
