function Data_ = auto_filter_LVP_cycles(data_,LV_outliers_threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LVP_ED_mode = 1;
RR_Interval_mode = 1;

dPdt_min_mode = 1;
dPdt_max_mode = 1;
 
Tau_mode = 1;

LVP_min_mode = 1;
LVP_max_mode = 1;

LV_outliers_mode = 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LVP_ED_percentiles =  [5,95];

RR_Interval_percentiles =  [5,95];


dPdt_min_percentiles =  [5,95];
dPdt_max_percentiles =  [5,95];

Tau_percentiles     =  [5,95];

LVP_min_percentiles =  [5,95];
LVP_max_percentiles =  [5,95];

%LV_outliers_threshold = 1.96;
%LV_outliers_threshold = 1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter data
dPdt_max_values  = cell2mat({data_.LVP_max}');
dPdt_min_values  = cell2mat({data_.LVP_min}');
LVP_ED_values   = cell2mat({data_.LVP_ED}');
RR_Interval_values = cell2mat({data_.RR_Interval}');
Tau_values = cell2mat({data_.Tau}');
LVP_min_values = cell2mat({data_.LVP_min}');
LVP_max_values = cell2mat({data_.LVP_max}');


LVP_ED_lower_band = prctile(LVP_ED_values,LVP_ED_percentiles(1));
LVP_ED_upper_band = prctile(LVP_ED_values,LVP_ED_percentiles(2));

RR_Interval_lower_band = prctile(RR_Interval_values,RR_Interval_percentiles(1));
RR_Interval_upper_band = prctile(RR_Interval_values,RR_Interval_percentiles(2));

dPdt_min_lower_band = prctile(dPdt_min_values ,dPdt_min_percentiles(1));
dPdt_min_upper_band = prctile(dPdt_min_values ,dPdt_min_percentiles(2));

dPdt_max_lower_band = prctile(dPdt_max_values ,dPdt_max_percentiles(1));
dPdt_max_upper_band = prctile(dPdt_max_values ,dPdt_max_percentiles(2));

Tau_lower_band = prctile(Tau_values ,Tau_percentiles(1));
Tau_upper_band = prctile(Tau_values ,Tau_percentiles(2));

LVP_min_lower_band = prctile(LVP_min_values ,LVP_min_percentiles(1));
LVP_min_upper_band = prctile(LVP_min_values ,LVP_min_percentiles(2));

LVP_max_lower_band = prctile(LVP_max_values ,LVP_max_percentiles(1));
LVP_max_upper_band = prctile(LVP_max_values ,LVP_max_percentiles(2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_1: LVP_ED range (10 percentile, median, 90 percentile)
if LVP_ED_mode==1
    for DD =1:size(LVP_ED_values,1)
        %if discretize(pR_flags(DD),[pR_interval_lower_band ,pR_interval_upper_band])==1 
        if (LVP_ED_lower_band<LVP_ED_values(DD)) && (LVP_ED_values(DD)<LVP_ED_upper_band)             
            LVP_ED_flags(DD) = 1;
        else
            LVP_ED_flags(DD) = 0;
        end
    end
else
    LVP_ED_flags = ones(1,size(data_,2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_2: RR_Interval range (10 percentile, median, 90 percentile)
if RR_Interval_mode==1
for EE =1:size(RR_Interval_values,1)
    %if discretize(pR_flags(DD),[pR_interval_lower_band ,pR_interval_upper_band])==1 
    if (RR_Interval_lower_band<RR_Interval_values(EE)) && (RR_Interval_values(EE)<RR_Interval_upper_band)             
        RR_Interval_flags(EE) = 1;
    else
        RR_Interval_flags(EE) = 0;
    end
end
else
    RR_Interval_flags = ones(1,size(data_,2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_3: dPdt_minrange (10 percentile, median, 90 percentile)
if dPdt_min_mode==1
for FF =1:size(dPdt_min_values,1)
    %if discretize(pR_flags(DD),[pR_interval_lower_band ,pR_interval_upper_band])==1 
    if (dPdt_min_lower_band<dPdt_min_values(FF)) && (dPdt_min_values(FF)<dPdt_min_upper_band)             
        dPdt_min_flags(FF) = 1;
    else
        dPdt_min_flags(FF) = 0;
    end
end
else
    dPdt_min_flags = ones(1,size(data_,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_4: dPdt_maz range (10 percentile, median, 90 percentile)
if dPdt_max_mode==1
for JJ =1:size(dPdt_max_values,1)
    %if discretize(pR_flags(DD),[pR_interval_lower_band ,pR_interval_upper_band])==1 
    if (dPdt_max_lower_band<dPdt_max_values(JJ)) && (dPdt_max_values(JJ)<dPdt_max_upper_band)             
        dPdt_max_flags(JJ) = 1;
    else
        dPdt_max_flags(JJ) = 0;
    end
end
else
   dPdt_max_flags = ones(1,size(data_,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_5: Tau range (10 percentile, median, 90 percentile)
if Tau_mode==1
for JJi =1:size(Tau_values,1)
    %if discretize(Tau_flags(DD),[Tau_lower_band ,Tau_upper_band])==1 
    if (Tau_lower_band<Tau_values(JJi)) && (Tau_values(JJi)<Tau_upper_band)             
        Tau_flags(JJi) = 1;
    else
        Tau_flags(JJi) = 0;
    end
end
else
     Tau_flags = ones(1,size(data_,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_6: LVP_min range (10 percentile, median, 90 percentile)
if LVP_min_mode==1
for Ji =1:size(LVP_min_values,1)
    %if discretize(LVP_min_flags(DD),[LVP_min_lower_band ,LVP_min_upper_band])==1 
    if (LVP_min_lower_band<LVP_min_values(Ji)) && (LVP_min_values(Ji)<LVP_min_upper_band)             
        LVP_min_flags(Ji) = 1;
    else
        LVP_min_flags(Ji) = 0;
    end
end
else
   LVP_min_flags = ones(1,size(data_,2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rule_7: LVP_max range (10 percentile, median, 90 percentile)
if LVP_max_mode==1
for iJ =1:size(LVP_max_values,1)
    %if discretize(LVP_max_flags(DD),[LVP_max_lower_band ,LVP_max_upper_band])==1 
    if (LVP_max_lower_band<LVP_max_values(iJ)) && (LVP_max_values(iJ)<LVP_max_upper_band)             
        LVP_max_flags(iJ) = 1;
    else
        LVP_max_flags(iJ) = 0;
    end
end
else
    LVP_max_flags = ones(1,size(data_,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exclude LV pressure traces outliers
outliers_LVP_flags = Exclude_LV_Pressures_Outliers(data_,LV_outliers_mode,LV_outliers_threshold);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Indxs = linspace(1,size(data_,2),size(data_,2))';
Flags = LVP_ED_flags'.*RR_Interval_flags'.*dPdt_min_flags'.*dPdt_max_flags'...
        .*Tau_flags'.*LVP_min_flags'.*LVP_max_flags'.*outliers_LVP_flags';

Data_ = data_;
zero_idxs = Flags.*Indxs;  zero_idxs = find(zero_idxs == 0);
Data_(zero_idxs)=[];




function outliers_LVP_flags = Exclude_LV_Pressures_Outliers(data,LV_outliers_mode,threshold)
    if LV_outliers_mode == 1
       %LVP_curves = {curve1; curve2; curve3; ...}; 
        curves = {data.LVP}';
        
        % Determine the maximum length among all curves
        maxLen = max(cellfun(@numel, curves));
        
        % Interpolate curves onto a common set of points
        interpCurves = NaN(length(curves), maxLen);
        for i = 1:length(curves)
            x = 1:numel(curves{i});
            xi = linspace(1, numel(curves{i}), maxLen);
            interpCurves(i, :) = interp1(x, curves{i}, xi, 'linear');
        end
        
        % Calculate average curve
        avgCurve = mean(interpCurves, 1);
        
        % Compute standard deviation for each point
        stdDev = std(interpCurves, 1);
        
        % Set a threshold for deviation
        %threshold = 1.96;  % Adjust this value based on your specific data
        
        % Identify outlier curves
        outliers = [];
        for i = 1:size(interpCurves, 1)
            deviations = abs(interpCurves(i, :) - avgCurve);
            if any(deviations > threshold * stdDev)
                outliers = [outliers; i];
            end
        end
        outliers_LVP_flags = ones(1,size(data,2));
        outliers_LVP_flags(outliers) = 0 ;
    else
        outliers_LVP_flags = ones(1,size(data,2));
    end
end

end

