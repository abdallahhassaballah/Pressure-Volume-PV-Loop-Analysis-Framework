function [Results_final,matched_volumes] = PV_loop_analysis (data_ed,cim,rr_match_threshold,LV_outliers_threshold,study_path)
L2 = 3;   L1 = 1.5;
FontSize = 26;
% Filter data
    data_ed = auto_filter_LVP_cycles(data_ed,LV_outliers_threshold);


    rr_volume = (60*1000)./cell2mat({cim.HR}');
    rr_pressure = cell2mat({data_ed.RR_Interval}');


    cim_match_summary = [];
    lineStyles = linspecer(size(cim,2));
    for jj=1:size(cim,2)
        cim_match_summary(jj).ID = cim(jj).ID;
        cim_match_summary(jj).model_name  = cim(jj).model_name;
        cim_match_summary(jj).ed_frame    = cim(jj).ed_frame;
        cim_match_summary(jj).es_frame    = cim(jj).es_frame;
        cim_match_summary(jj).nframes     = cim(jj).nframes;
        cim_match_summary(jj).EDV         = cim(jj).EDV;
        cim_match_summary(jj).ESV         = cim(jj).ESV; 
        cim_match_summary(jj).SV          = cim(jj).SV;
        cim_match_summary(jj).MASS        = cim(jj).MASS;
        cim_match_summary(jj).EF          = cim(jj).EF;
        cim_match_summary(jj).GLS         = cim(jj).GLS;
        cim_match_summary(jj).HR          = cim(jj).HR;
        cim_match_summary(jj).Volumes     = cim(jj).Volumes;
        cim_match_summary(jj).rr_pressure_match_idx = find_matcing_rr_pressure_cycles(rr_pressure,rr_volume(jj),rr_match_threshold);
        cim_match_summary(jj).cim_color = lineStyles(jj,:);
    end



    cim_idxs = {cim_match_summary.rr_pressure_match_idx}';
    emptyRows = find(cellfun(@isempty, cim_idxs));

    cim_match_summary_new = cim_match_summary;
    cim_match_summary_new(emptyRows) = [];

    matched_volumes = cim_match_summary_new;
    %% Plot CIM Volume Curves
    figure; hold on
    for jj=1:size(cim,2)
        Volume = upsampling_CIM_volume(cim(jj).Volumes',200);
        time_vol  = linspace(0,60/cim(jj).HR,size(Volume,1))'*1000;
        if ismember(jj, emptyRows)==1
            plot(time_vol,Volume,'--','color',[0.5 0.5 0.5],'LineWidth',L2);
        else
            plot(time_vol,Volume,'-','color',lineStyles(jj,:),'LineWidth',L2);
        end
    end
    xlabel('Time (ms)','FontWeight','bold');
    ylabel('Volume (mL)','FontWeight','bold');
    box on; set(gcf,'color','w');
    set(gca,'FontSize',FontSize, 'FontName','Times New Roman')
    set(gcf, 'PaperUnits', 'inches');
    cd(study_path)
    set(gcf, 'PaperPosition', [0 0 12 8]);
    print(gcf,'CIM_Volume_Curves.png','-dpng'); 

    %% Plot_Selected_ECG_LVP_Traces
    Plot_Selected_ECG_LVP_Traces(data_ed)
    cd(study_path)
    set(gcf, 'PaperPosition', [0 0 12 8]);
    print(gcf,'Selected_ECG_LVP_Traces.png','-dpng'); 

     %% Plot LV PV loops
    smooth_factor = 2;
    outlier_factor = 5;
    
    Results_final = [];
   
    fig = figure; hold on;
    for jk =1:size(cim_match_summary_new,2)   
    cim_summary = cim_match_summary_new(jk);
    
    cim_sumary_ed = Extract_cim_Summary(data_ed(cim_summary.rr_pressure_match_idx),cim_summary);
    
    idxs = cell2mat({cim_summary.rr_pressure_match_idx}');
    data_ed__= data_ed(idxs);
        for kk = 1:size(idxs,1)       
            LVP_avg_ = upsampling_CIM_volume(data_ed__(kk).LVP,200);
            vol_avg_ = upsampling_CIM_volume(cim_summary.Volumes',200);
            t = linspace(0,length(LVP_avg_)/240,length(LVP_avg_));
            cim_color = cim_summary.cim_color;
            [p_no_outliers,v_no_outliers,t_no_outliers] = plot_pv_loop_smoothed(LVP_avg_', vol_avg_', t, smooth_factor, outlier_factor,cim_color);
            ESV = cim_summary.ESV;
            SV = cim_summary.SV;
            HR = data_ed__(kk).HR; 
            [SW, PE, PVA, eta,EEV, CO, Ea,Ees,Coupling_ratio,ESP,MEP] = Extract_PV_Loop_Paramaters(p_no_outliers,v_no_outliers,ESV,SV,HR);
        
            data_ed__(kk).cim_model_name = cim_summary.model_name;
            data_ed__(kk).EDV = cim_summary.EDV; 
            data_ed__(kk).ESV = cim_summary.ESV; 
            data_ed__(kk).SV = cim_summary.SV;
            data_ed__(kk).LVMass = cim_summary.MASS;
            data_ed__(kk).EF = cim_summary.EF;
            data_ed__(kk).GLS = cim_summary.GLS;
            data_ed__(kk).cim_HR = cim_summary.HR;
            data_ed__(kk).cim_ed_frame = cim_summary.ed_frame; 
            data_ed__(kk).cim_es_frame = cim_summary.es_frame;
            data_ed__(kk).cim_nframes = cim_summary.nframes;
            data_ed__(kk).cim_volumes = cim_summary.Volumes;
            data_ed__(kk).cim_color = cim_summary.cim_color;
            
            data_ed__(kk).p_filtered = p_no_outliers;
            data_ed__(kk).v_filtered = v_no_outliers;
            data_ed__(kk).t_filtered = t_no_outliers;
            data_ed__(kk).SW =SW;
            data_ed__(kk).PE=PE;
            data_ed__(kk).PVA=PVA;
            data_ed__(kk).eta=eta;
            data_ed__(kk).EEV=EEV;
            data_ed__(kk).CO=CO;
            data_ed__(kk).Ea=Ea;
            data_ed__(kk).Ees=Ees;
            data_ed__(kk).Coupling_ratio=Coupling_ratio;
            data_ed__(kk).ESP=ESP;
            data_ed__(kk).MEP=MEP;
        end
      Results_final = [Results_final, data_ed__];
    end
    xlim auto
    ylim auto
    cd(study_path)
    set(gcf, 'PaperPosition', [0 0 12 8]);
    print(gcf,'PV_Loops.png','-dpng'); 
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Internal Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = find_matcing_rr_pressure_cycles(rr_pressure_values,rr_volume_match_value, threshold_percentage)

    threshold_value = rr_volume_match_value*(threshold_percentage/100);
    min_value = rr_volume_match_value - threshold_value;
    max_value = rr_volume_match_value + threshold_value;

    
    indices = find(rr_pressure_values >= min_value & rr_pressure_values <= max_value);

end



function [p_no_outliers,v_no_outliers,t_no_outliers] = plot_pv_loop_smoothed(p, v, t,smooth_factor, outlier_factor,cim_color)
L2 = 3;   L1 = 1.5;  sz = 750;
FontSize = 26;
color = [0.904705882352941,0.191764705882353,0.198823529411765];
p_ =p;
v_ =v;

  % Ensure that the end-diastolic points have the same values
    p(1) = (p_(1)+p_(end))./2;
    v(1) = (v_(1)+v_(end))./2;

    p(end) = (p_(1)+p_(end))./2;
    v(end) = (v_(1)+v_(end))./2;
    

% p: vector matrix of pressure values
% v: vector matrix of volume values
% smooth_factor: smoothing factor (must be an odd integer)
% outlier_factor: factor for removing outliers (values outside of the range of median +/- outlier_factor*IQR)

% Check that p and v are the same size
if ~isequal(size(p), size(v))
    error('p and v must be the same size.')
end

% Check that p and v are vectors
if ~(isvector(p) && isvector(v))
    error('p and v must be vectors.')
end

% Ensure that the PV loop is closed and starts and ends at the same point
if p(1) ~= p(end) || v(1) ~= v(end)
    p = [p, p(1)];
    v = [v, v(1)];
end

% Check that the resulting vectors are the same size
if ~isequal(size(p), size(v))
    error('Error closing the PV loop.')
end

% Apply smoothing to the pressure and volume data
p_smooth = smoothdata(p, 'sgolay', smooth_factor);
v_smooth = smoothdata(v, 'sgolay', smooth_factor);

% Remove outliers from the smoothed data
p_median = median(p_smooth);
p_iqr = iqr(p_smooth);
v_median = median(v_smooth);
v_iqr = iqr(v_smooth);

p_min = p_median - outlier_factor*p_iqr;
p_max = p_median + outlier_factor*p_iqr;
v_min = v_median - outlier_factor*v_iqr;
v_max = v_median + outlier_factor*v_iqr;

p_outliers = p_smooth < p_min | p_smooth > p_max;
v_outliers = v_smooth < v_min | v_smooth > v_max;
outliers = p_outliers | v_outliers;

p_no_outliers = p_smooth(~outliers);
v_no_outliers = v_smooth(~outliers);

% Plot the PV loop with smoothed and outlier-removed data
%scatter(v_no_outliers, p_no_outliers,sz,'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerFaceAlpha',0.75);
%plot(v_, p_, '-k','LineWidth',L2);
plot(v_no_outliers, p_no_outliers, '-','color',cim_color,'LineWidth',L2);
hold on;
xlabel('Volume (mL)','FontWeight','bold');
ylabel('Pressure (mmHg)','FontWeight','bold');
%title('Smoothed and Outlier-Removed PV Loop')
set(gcf,'color','w');
set(gca,'FontSize',FontSize, 'FontName','Times New Roman')
set(gcf, 'PaperUnits', 'inches');
box on

% Set axis limits
xlim([min(v)-0.1*(max(v)-min(v)), max(v)+0.1*(max(v)-min(v))])
ylim([min(p)-0.1*(max(p)-min(p)), max(p)+0.1*(max(p)-min(p))])

% Ensure that the PV loop appears as a reasonable rectangular shape
ratio = (max(p_no_outliers)-min(p_no_outliers))/(max(v_no_outliers)-min(v_no_outliers));
if ratio > 1.5 || ratio < 0.5
    warning('The PV loop may not appear as a reasonable rectangular shape.')
end


t_no_outliers = linspace(t(1),t(2), length(p_no_outliers));

end


function Plot_Selected_ECG_LVP_Traces(data)

addpath(genpath('\\nas1\hpc\Heart_Biomechanics_Files\MATLAB_Packages'));
C = linspecer(10,'qualitative');
C_ = linspecer(11,'qualitative');
L2 = 3;   L1 = 1.5;
FontSize = 26;
    figure; hold on
    for jj=1:size(data,2)
        yyaxis left;
        plot(data(jj).time*1000,data(jj).LVP,'-','color','k','LineWidth',L2);
        scatter_new(data(jj).time*1000,data(jj).LVP,data(jj).ED_idx,C(5,:));
        scatter_new(data(jj).time*1000,data(jj).LVP,data(jj).DS_idx,C_(11,:));
        yyaxis right;
        plot(data(jj).time*1000,data(jj).ECG,':','color', [0.904705882352941,0.191764705882353,0.198823529411765],'linewidth', L1);
        %{
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).R_Peak);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).P_Peak);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).T_Peak);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).S_Peak);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).Q_Peak);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).P_Offset);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).T_Offset);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).P_Onset);
        scatter_new(data(jj).time*1000,data(jj).ECG,data(jj).T_Onset);
        %}
    end
    xlabel('Time (ms)','FontWeight','bold');
    yyaxis left;    ylabel('Pressure (mmHg)','FontWeight','bold');
    yyaxis right;    ylabel('ECG (mV)','FontWeight','bold');
    box on; set(gcf,'color','w');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = [0.904705882352941,0.191764705882353,0.198823529411765];
    set(gca,'FontSize',FontSize, 'FontName','Times New Roman')
    set(gcf, 'PaperUnits', 'inches');

    function scatter_new(time,ECG,idx,color)
        sz = 750;
        if ~isnan(idx)==1
            scatter(time(idx),ECG(idx),sz,'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerFaceAlpha',0.75);
        end
    end
end



function integrate_PV_area(PV)

% Load the data
% PV = load('your_data.mat');  % Uncomment this if your data is in a .mat file

% Assume the first row is volume (V) and second row is pressure (P)
V = PV(1,:);
P = PV(2,:);

% Ensure that the loop is closed. If it is not, close it.
if V(1) ~= V(end) || P(1) ~= P(end)
    V = [V, V(1)];
    P = [P, P(1)];
end

% Apply the trapezoidal rule for numerical integration
Area = trapz(V, P);

% Display the result
disp(['The area under the PV loop is: ', num2str(Area), ' units^2'])

end

end