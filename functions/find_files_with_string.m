function file_paths = find_files_with_string(main_dir, search_str)

% Get a list of all files and folders in the main directory
all_files = dir(main_dir);

% Initialize a cell array to store the file paths
file_paths = {};

% Loop through all the files and folders in the main directory
for i = 1:length(all_files)
    file = all_files(i);
    
    % Check if the current item is a directory
    if file.isdir
        % Skip the '.' and '..' directories
        if strcmp(file.name, '.') || strcmp(file.name, '..')
            continue
        end
        
        % Recursively search for files in the subdirectory
        subdir_path = fullfile(main_dir, file.name);
        subdir_files = find_files_with_string(subdir_path, search_str);
        
        % Add the subdirectory files to the list of file paths
        file_paths = [file_paths, subdir_files];
    else
        % Check if the current file name contains the search string
        if ~isempty(strfind(file.name, search_str))
            % If the file name contains the search string, add its path to the list
            file_paths{end+1} = fullfile(main_dir, file.name);
        end
    end
end
