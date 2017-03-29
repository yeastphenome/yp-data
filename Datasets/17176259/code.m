%% Bishop~Avery, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bishop_avery_2007.pmid = 17176259;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bishop_avery_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mata_DeletionArray+slow_growers.xlsx', '96');
tested_orfs = tested.raw(3:end,2);
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mata_DeletionArray+slow_growers.xlsx', 'slow growers');
tested_orfs = [tested_orfs; tested.raw(3:end,2)];

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_genenames_resistant] = read_data('textread','./raw_data/hits_genenames_resistant.txt', '%s');

hits_genenames_resistant = clean_genename(hits_genenames_resistant);
hits_orfs_resistant = translate(hits_genenames_resistant);
hits_orfs_resistant = unique(hits_orfs_resistant);

hits_scores_resistant = ones(length(hits_orfs_resistant),1);

[FILENAMES{end+1}, hits_genenames_sensitive] = read_data('textread','./raw_data/hits_genenames_sensitive.txt', '%s');

hits_genenames_sensitive = clean_genename(hits_genenames_sensitive);
hits_orfs_sensitive = translate(hits_genenames_sensitive);
hits_orfs_sensitive = unique(hits_orfs_sensitive);

hits_scores_sensitive = -ones(length(hits_orfs_sensitive),1);

overlapping = intersect(hits_orfs_resistant, hits_orfs_sensitive);

hits_orfs = [hits_orfs_resistant; hits_orfs_sensitive];
hits_scores = [hits_scores_resistant; hits_scores_sensitive];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 5 orfs to be added

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [173];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
bishop_avery_2007.orfs = tested_orfs;
bishop_avery_2007.ph = hit_data_names;
bishop_avery_2007.data = zeros(length(bishop_avery_2007.orfs),length(bishop_avery_2007.ph));
bishop_avery_2007.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, bishop_avery_2007.orfs);
bishop_avery_2007.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./bishop_avery_2007.mat','bishop_avery_2007');

%% Print out

fid = fopen('./bishop_avery_2007.txt','w');
write_matrix_file(fid, bishop_avery_2007.orfs, bishop_avery_2007.ph, bishop_avery_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bishop_avery_2007)
end

end
