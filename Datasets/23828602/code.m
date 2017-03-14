%% Ding~Bakalinsky, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ding_bakalinsky_2013.pmid = 23828602;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ding_bakalinsky_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES, hits_genenames_R] = read_data('textscan', './raw_data/hits_genenames_R.txt', '%s');

% Clean the names
hits_genenames_R = clean_genename(hits_genenames_R);

% Translate from genenames to ORFs
hits_orfs_R = translate(hits_genenames_R);

% Get the unique length
hits_orfs_R = unique(hits_orfs_R);

% Make the hit strains data
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

% Make sure the that all the hits are part of the tested set
[missing, ix] = setdiff(hits_orfs_R, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [144];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
ding_bakalinsky_2013.orfs = tested_orfs;
ding_bakalinsky_2013.ph = hit_data_names;
ding_bakalinsky_2013.data = zeros(length(ding_bakalinsky_2013.orfs),length(ding_bakalinsky_2013.ph));
ding_bakalinsky_2013.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs_R, ding_bakalinsky_2013.orfs);
ding_bakalinsky_2013.data(ind2) = hits_scores_R(ind1);

%% Save

save('./ding_bakalinsky_2013.mat','ding_bakalinsky_2013');

%% Print out

fid = fopen('./ding_bakalinsky_2013.txt','w');
write_matrix_file(fid, ding_bakalinsky_2013.orfs, ding_bakalinsky_2013.ph, ding_bakalinsky_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ding_bakalinsky_2013)
end

end