%% Auesukaree~Harashima, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
auesukaree_harashima_2009.pmid = 19638689;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(auesukaree_harashima_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

treatments = {'EtOH';'MeOH';'propanol';'NaCl';'H2O2';'37C'};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat alpha_KOset list.xlsx');
tested_orfs = tested.raw(4:end,2);

tested_orfs = clean_orf(tested_orfs);
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Load data

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [162 432 433 434 435 436]';

[FILENAMES{end+1}, data_hits{1}] = read_data('textread','./raw_data/ethanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{2}] = read_data('textread','./raw_data/methanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{3}] = read_data('textread','./raw_data/propanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{4}] = read_data('textread','./raw_data/nacl_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{5}] = read_data('textread','./raw_data/h2o2_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{6}] = read_data('textread','./raw_data/heat_sensitivity_hits.txt', '%s');

data_hits_orfs = cell(size(data_hits));

data_hits_orfs{1} = unique(translate(data_hits{1}));
[missing, ix] = setdiff(data_hits_orfs{1}, tested_orfs);
data_hits_orfs{1}(strcmp('YHR039C-A', data_hits_orfs{1})) = {'YHR039C-B'};

data_hits_orfs{2} = unique(translate(data_hits{2}));
[missing, ix] = setdiff(data_hits_orfs{2}, tested_orfs);

data_hits_orfs{3} = unique(translate(data_hits{3}));
[missing, ix] = setdiff(data_hits_orfs{3}, tested_orfs);

data_hits_orfs{4} = unique(translate(data_hits{4}));
[missing, ix] = setdiff(data_hits_orfs{4}, tested_orfs);

data_hits_orfs{5} = unique(translate(data_hits{5}));
[missing, ix] = setdiff(data_hits_orfs{5}, tested_orfs);

data_hits_orfs{6} = unique(translate(data_hits{6}));
[missing, ix] = setdiff(data_hits_orfs{6}, tested_orfs);
tested_orfs = [tested_orfs; missing];

hit_data = zeros(length(tested_orfs), length(hit_data_ids));
for i = 1 : length(hit_data_ids)
    [~,~,ind2] = intersect(data_hits_orfs{i}, tested_orfs);
    hit_data(ind2,i) = -1;
end

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

auesukaree_harashima_2009.orfs = tested_orfs;
auesukaree_harashima_2009.data = hit_data;
auesukaree_harashima_2009.ph = hit_data_names;
auesukaree_harashima_2009.dataset_ids = hit_data_ids;

%% Save

save('./auesukaree_harashima_2009.mat','auesukaree_harashima_2009');

%% Print out

fid = fopen('./auesukaree_harashima_2009.txt','w');
write_matrix_file(fid, auesukaree_harashima_2009.orfs, auesukaree_harashima_2009.ph, auesukaree_harashima_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(auesukaree_harashima_2009)
end

end
