%% Gardocki~Lopes, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gardocki_lopes_2005.pmid = 15755922;

phenotypes = {'expression of PIS1'};
treatments = {'glucose 2%';'glycerol 3%'};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/120_GENES_AFFECTING_PIS1.xlsx');

orfs = data.raw(8:end,2);
raw_data = data.raw(8:end,[11 14]);

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];
raw_data(inds,:) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
raw_data(inds,:) = [];

% Transform symbols to numbers
% Glucose: + -> -1, - -> 0
% Glycerol: - -> 0, + -> +1

inds = find(strcmp('+', raw_data(:,1)));
raw_data(inds,1) = {-1};
inds = find(strcmp('-', raw_data(:,1)));
raw_data(inds,1) = {0};
inds = find(strcmp('+', raw_data(:,2)));
raw_data(inds,2) = {1};
inds = find(strcmp('-', raw_data(:,2)));
raw_data(inds,2) = {0};
inds = find(strcmp('NT', raw_data));
raw_data(inds) = {NaN};

raw_data = cell2mat(raw_data);

orfs = upper(strtrim(orfs));
[t,t2] = grpstats(raw_data, orfs,{'mean','gname'});

gardocki_lopes_2005.orfs = t2;
gardocki_lopes_2005.data = t;
gardocki_lopes_2005.ph = strcat(phenotypes, '; ', treatments);

save('./gardocki_lopes_2005.mat','gardocki_lopes_2005');

fid = fopen('./gardocki_lopes_2005.txt','w');
write_matrix_file(fid, gardocki_lopes_2005.orfs, gardocki_lopes_2005.ph, gardocki_lopes_2005.data);
fclose(fid);

end
