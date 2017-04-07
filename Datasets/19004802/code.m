%% Ralser~Lehrach, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ralser_lehrach_2008.pmid = 19004802;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ralser_lehrach_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES{end+1}, hits_orfs] = read_data('textread','./raw_data/hits_orfs.txt', '%s');

hits_orfs = clean_orf(hits_orfs);

% If in gene name form, transform into ORF name
hits_orfs = translate(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));  

hits_orfs = unique(hits_orfs);
hits_data = ones(size(hits_orfs));

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v2.0.xlsx');
tested_orfs = tested.raw(2:end,1);

tested_orfs = clean_orf(tested_orfs);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs  = translate(tested_orfs);

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs = unique(tested_orfs);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [552];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
ralser_lehrach_2008.orfs = tested_orfs;
ralser_lehrach_2008.ph = hit_data_names;
ralser_lehrach_2008.data = zeros(length(ralser_lehrach_2008.orfs),length(ralser_lehrach_2008.ph));
ralser_lehrach_2008.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, ralser_lehrach_2008.orfs);
ralser_lehrach_2008.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./ralser_lehrach_2008.mat','ralser_lehrach_2008');

%% Print out

fid = fopen('./ralser_lehrach_2008.txt','w');
write_matrix_file(fid, ralser_lehrach_2008.orfs, ralser_lehrach_2008.ph, ralser_lehrach_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ralser_lehrach_2008)
end

end
