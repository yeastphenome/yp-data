%% Copic~Miller, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
copic_miller_2009.pmid = 19433630;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(copic_miller_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/mat_alpha_obs_v1.0.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,1);

tested_orfs = clean_orf(tested_orfs);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

tested_orfs = translate(tested_orfs);

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds))

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS1.xlsx', 'TABLE S1');
hits_orfs = data.raw(5:end,1);
hits_notes = data.raw(5:end,3);
hits_scores = data.raw(5:end,6);

inds = strcmp('mat-a only', hits_notes);
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = clean_orf(hits_orfs);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = translate(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds))

hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_scores = cell2mat(hits_scores);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [194];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
copic_miller_2009.orfs = tested_orfs;
copic_miller_2009.ph = hit_data_names;
copic_miller_2009.data = zeros(length(copic_miller_2009.orfs),length(copic_miller_2009.ph));
copic_miller_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, copic_miller_2009.orfs);
copic_miller_2009.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./copic_miller_2009.mat','copic_miller_2009');

%% Print out

fid = fopen('./copic_miller_2009.txt','w');
write_matrix_file(fid, copic_miller_2009.orfs, copic_miller_2009.ph, copic_miller_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(copic_miller_2009)
end

end
