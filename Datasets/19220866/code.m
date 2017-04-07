%% Mira~Sa-Correia, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mira_sa_correia_2009.pmid = 19220866;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mira_sa_correia_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested

[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

tested_orfs = clean_orf(tested_orfs);

tested_orfs = translate(tested_orfs);

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds))

tested_orfs = unique(tested_orfs);

%% Load data

[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');
[FILENAMES{end+1}, hits_genenames_moderate] = read_data('textread','./raw_data/hits_genenames_moderate.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_genenames_moderate = clean_genename(hits_genenames_moderate);

hits_genenames_strong = setdiff(hits_genenames, hits_genenames_moderate);

hits_orfs_moderate = translate(hits_genenames_moderate);
hits_orfs_moderate = unique(hits_orfs_moderate);

[hits_orfs_strong, translated] = translate(hits_genenames_strong);
hits_orfs_strong(~translated & ~is_orf(hits_orfs_strong)) = [];
hits_orfs_strong = unique(hits_orfs_strong);

[missing, ix] = setdiff(hits_orfs_strong, tested_orfs);
tested_orfs = [tested_orfs; missing];      % 3 ORFs added to tested

[missing, ix] = setdiff(hits_orfs_moderate, tested_orfs);

hits_data_strong = zeros(length(hits_orfs_strong),1)-2;
hits_data_moderate = zeros(length(hits_orfs_moderate),1)-1;

hits_orfs = [hits_orfs_strong; hits_orfs_moderate];
hits_data = [hits_data_strong; hits_data_moderate];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [157];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
mira_sa_correia_2009.orfs = tested_orfs;
mira_sa_correia_2009.ph = hit_data_names;
mira_sa_correia_2009.data = zeros(length(mira_sa_correia_2009.orfs),length(mira_sa_correia_2009.ph));
mira_sa_correia_2009.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, mira_sa_correia_2009.orfs);
mira_sa_correia_2009.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./mira_sa_correia_2009.mat','mira_sa_correia_2009');

%% Print out

fid = fopen('./mira_sa_correia_2009.txt','w');
write_matrix_file(fid, mira_sa_correia_2009.orfs, mira_sa_correia_2009.ph, mira_sa_correia_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mira_sa_correia_2009)
end

end
