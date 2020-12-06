%% Gardocki~Lopes, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gardocki_lopes_2005.pmid = 15755922;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gardocki_lopes_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/120_GENES_AFFECTING_PIS1.xlsx');

orfs = data.raw(8:end,2);
raw_data = data.raw(8:end,[11 14]);

orfs = clean_orf(orfs);

orfs(find(strcmp('YYKL138C', orfs))) = {'YKL138C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(orfs));
disp(orfs(inds)); 

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

[raw_data,orfs] = grpstats(raw_data, orfs,{'mean','gname'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [575 576]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
gardocki_lopes_2005.orfs = orfs;
gardocki_lopes_2005.ph = hit_data_names;
gardocki_lopes_2005.data = raw_data;
gardocki_lopes_2005.dataset_ids = hit_data_ids;

%% Save

save('./gardocki_lopes_2005.mat','gardocki_lopes_2005');

%% Print out

fid = fopen('./gardocki_lopes_2005.txt','w');
write_matrix_file(fid, gardocki_lopes_2005.orfs, gardocki_lopes_2005.ph, gardocki_lopes_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gardocki_lopes_2005)
end

end
