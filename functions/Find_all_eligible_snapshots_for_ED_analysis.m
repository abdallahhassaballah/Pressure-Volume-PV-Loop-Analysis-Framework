function [file_paths_ECG,file_paths_LV]= Find_all_eligible_snapshots_for_ED_analysis(main_dir)


search_str = '_h5.h5';
file_paths = find_files_with_string(main_dir, search_str);
file_paths_ECG = file_paths';


% Define the old and new strings
old_string = '_h5.h5';
new_string = '.mat';

% Replace the string in the cell array
cell_array_ = strrep(file_paths, old_string, new_string);

% Define the old and new strings
old_string = '\ECG_LV_results';
new_string = '\Data';

% Replace the string in the cell array
cell_array = strrep(cell_array_, old_string, new_string);


file_paths_LV = cell_array';


end


