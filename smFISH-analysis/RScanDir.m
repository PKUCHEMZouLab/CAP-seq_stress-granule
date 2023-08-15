% Get the all file with specific string (file_mask, like '*.tif') in selected path (path)
% to a cell variants named file_list
% Zoulab, 2020
function [file_list] = RScanDir(path, file_mask)
p = genpath(path);      % get all subfiles in the path and write into variant p, splited by a ';'
length_p = size(p,2);   % get the length of p
path = {};
temp = [];
file_list = {};
path_num = 1;
for i = 1:length_p %Search for the delimiter ';' and once found, write the path temp into the path_temp array.
    if p(i) ~= ';'
        temp = [temp p(i)];
    else 
        temp = [temp '\']; % add '\' at the end of the path
        path = [path ; temp];
        temp = [];
    end
end  
clear p length_p temp;
% pickup the 
file_num = size(path,1);% the number of subfolders
for i = 1:file_num
    file_path =  path{i};
    img_path_list = dir(strcat(file_path, file_mask)); % search the files with specific string in each subfolder
    img_num = length(img_path_list); % get the number of selected files in current subfolder
    if img_num > 0
        for j = 1:img_num
            image_name = img_path_list(j).name; % get the filename with extension.
%             fprintf('We get %s\n', strcat(file_path,image_name));    % input the real path of selected file
            file_list{path_num} = [strcat(file_path,image_name)];
            path_num = path_num+1
        end
    end
end
end
