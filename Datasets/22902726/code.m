%% Schmidt~Boyer, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
schmidt_boyer_2012.pmid = 22902726;

phenotypes = {'growth (OD)';'growth (colony size)';'growth (MIC)'};
treatments = {'boric acid'};

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(schmidt_boyer_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/BA sensitivity comprehensive data file.xlsx', 'Strain list');
tested_orfs = tested.raw(2:end,2);

tested_orfs = clean_orf(tested_orfs);

inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/BA sensitivity comprehensive data file.xlsx', 'Sheet1');

hit_orfs = data.raw(3:end,1);

hit_orfs = clean_orf(hit_orfs);

inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];

hit_orfs = unique(hit_orfs);
hit_data = -ones(length(hit_orfs),1);

[missing, ix] = setdiff(hit_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [424];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

schmidt_boyer_2012.orfs = tested_orfs;
schmidt_boyer_2012.data = zeros(length(tested_orfs), length(hit_data_names));
[~,ind1,ind2] = intersect(hit_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,1) = hit_data(ind1);
schmidt_boyer_2012.ph = hit_data_names;
schmidt_boyer_2012.dataset_ids = hit_data_ids;

%% Save

save('./schmidt_boyer_2012.mat','schmidt_boyer_2012');


%% Print out

fid = fopen('./schmidt_boyer_2012.txt','w');
write_matrix_file(fid, schmidt_boyer_2012.orfs, schmidt_boyer_2012.ph, schmidt_boyer_2012.data);
fclose(fid);

end
