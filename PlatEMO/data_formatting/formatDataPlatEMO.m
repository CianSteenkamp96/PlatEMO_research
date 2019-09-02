% Cian Steenkamp
% This script creates files where each line represents a found solution; columns seperated by spaces (column for each objective).
% This script formats the 'archived' solutions found during the search.
% This script will need to be run in each
% opt_raw_data/algo/instance/problem/objective/file (for KnEA, MOEADD, CDAS_SMPSO) folder
% and
% opt_raw_data/algo/problem/objective/file (for NSGAIII) folder.
% Then move the created files to the appropriate subfolders in the:
% opt_formatted_data/ folder.
								
clc;
clear;

files = dir('*.mat');

for file = files'
    d = load(file.name);
    s = size(d.result, 1);
    for i = (1 + s): (2 * s)
        disp(d.result{i}.objs);
        splitted = split(file.name, '.');
        fileName = string(strcat(splitted(1, 1), '_iter_', int2str(d.result{i - s})));
        dlmwrite(fileName, d.result{i}.objs, 'precision', '%.4f', 'delimiter', ' ');
    end
end
