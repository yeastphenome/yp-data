%% Cheng~Bakalinsky, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
cheng_bakalinsky_2007.pmid = 17644632;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(cheng_bakalinsky_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
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

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supplementary_Table_1.xlsx', 'Table 1');

% Dataset1: Tested = all; hits = liquid
hits_orfs = data.raw(4:end,2);
hits_scores = data.raw(4:end,3);

hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));  

hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_scores = cellfun(@str2num, hits_scores);

% If the same strain is present more than once, average its values
[hits_orfs, hits_scores] = grpstats(hits_scores, hits_orfs, {'gname','mean'});

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [114];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
cheng_bakalinsky_2007.orfs = tested_orfs;
cheng_bakalinsky_2007.ph = hit_data_names;
cheng_bakalinsky_2007.data = zeros(length(cheng_bakalinsky_2007.orfs),length(cheng_bakalinsky_2007.ph));
cheng_bakalinsky_2007.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, cheng_bakalinsky_2007.orfs);
cheng_bakalinsky_2007.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./cheng_bakalinsky_2007.mat','cheng_bakalinsky_2007');

%% Print out

fid = fopen('./cheng_bakalinsky_2007.txt','w');
write_matrix_file(fid, cheng_bakalinsky_2007.orfs, cheng_bakalinsky_2007.ph, cheng_bakalinsky_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(cheng_bakalinsky_2007)
end

end
