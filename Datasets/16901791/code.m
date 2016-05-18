%% Parsons~Boone, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
parsons_boone_2006.pmid = 16901791;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(parsons_boone_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load hit strains and tested strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/mmc8.xlsx', 'Supplementary Table 7');

% Get list of ORFs from both sets of ORFs
hit_orfs = data(4:end,1);

% Get Treatments
treatments = strtrim(data(3, 2:end));

% Get rid of empty rows
inds = find(cellfun(@isnumeric, hit_orfs));
hit_orfs(inds) = [];

% Clean up ORFs
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds)); 

% Retrieve data
hit_data = data(4:end, 2:end);

inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Data Calculation.
hit_data = -hit_data;

% Average any repeated value
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[FILENAMES{end+1}, dataset_map] = read_data('readtable', './extras/phenotype_dataset.txt', 'Delimiter','\t', 'ReadVariableNames', false);

[~,ind1,ind2] = intersect(treatments, strtrim(dataset_map.Var1));
hit_data_ids = nan(size(treatments));
hit_data_ids(ind1) = dataset_map.Var2(ind2);

[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

parsons_boone_2006.orfs = hit_orfs;
parsons_boone_2006.ph = hit_data_names;
parsons_boone_2006.data = hit_data;
parsons_boone_2006.dataset_ids = hit_data_ids; 

%% Save

save('./parsons_boone_2006.mat','parsons_boone_2006');

fid = fopen('./parsons_boone_2006.txt','w');
write_matrix_file(fid, parsons_boone_2006.orfs, parsons_boone_2006.ph, parsons_boone_2006.data);
fclose(fid);

end
