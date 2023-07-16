function number_of_LV_files = load_invasive_raw_pressure_ecg_data(pressure_path,result_path)
    cd(pressure_path)
    List = dir('*.inf');
    N = 1;
    for i = 1:size(List,1)
        fid=fopen(List(i).name);
        file = regexp(fileread(List(i).name),'\n','split');
        if contains(file{1,2},'LV')
            file(1:10) = []; 
            LeadI_idx = find(contains(file,'II'));  LeadI_idx = LeadI_idx(1);
            LV_idx = find(contains(file,'LV'));    LV_idx    = LV_idx(end);

            T = readtable(append(erase(List(i).name,'.inf'),'.txt'));
   
            Data = [smooth(table2array(T(:,2))),smooth(table2array(T(:,LV_idx)))];
            ECG = array2table(smooth(table2array(T(:,2)))); ECG.Properties.VariableNames(1) = {'ECG'};
            LVP = array2table(smooth(table2array(T(:,LV_idx))));  LVP.Properties.VariableNames(1) = {'LVP'};

            if ~exist(result_path, 'dir')
                mkdir(result_path)
            end

            cd(result_path)
            writetable(ECG,append('ECG_LV',num2str(N),'.txt'),'Delimiter',',')
            save(append('Data_',num2str(N),'.mat'),'Data')
            N = N+1;

            cd(pressure_path)       
        end
    end
 cd(result_path)
 number_of_LV_files = N-1;
 Summary.number_of_LV_files= N-1;
 T = struct2table(Summary);
 writetable(T,'Summary.csv','Delimiter',',')
 type 'Summary.csv';
end