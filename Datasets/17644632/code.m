%% Cheng~Bakalinsky, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

cheng_bakalinsky_2007.pmid = 17644632;

phenotypes = {'growth (MIC)'};
treatments = {'oxalic acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/YSC1054Y.copy.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supplementary_Table_1.xlsx', 'Table 1');

% Dataset1: Tested = all; hits = liquid
hits_orfs = data.raw(4:end,2);
hits_scores = data.raw(4:end,3);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cellfun(@str2num, hits_scores);

% Many orfs are listed multiple times, but their scores are always the
% same. So, I can easily take the first of their values

[hits_orfs, ia,ic] = unique(hits_orfs);
hits_scores = hits_scores(ia);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

cheng_bakalinsky_2007.orfs = tested_orfs;
cheng_bakalinsky_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
cheng_bakalinsky_2007.data(ind2,1) = hits_scores(ind1);

cheng_bakalinsky_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('./cheng_bakalinsky_2007.mat','cheng_bakalinsky_2007');
return;

% Save data into database
dt = cheng_bakalinsky_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./cheng_bakalinsky_2007.txt','w');
write_matrix_file(fid, cheng_bakalinsky_2007.orfs, cheng_bakalinsky_2007.ph, cheng_bakalinsky_2007.data);
fclose(fid);

end

