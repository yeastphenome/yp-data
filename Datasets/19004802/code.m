%% Ralser~Lehrach, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ralser_lehrach_2008.pmid = 19004802;

phenotypes = {'growth'};
treatments = {'2-DG'};

% Load data
[FILENAMES{end+1}, hits_orfs] = dataread('textread','./raw_data/hits_orfs.txt', '%s');
hits_orfs = unique(strtrim(upper(hits_orfs)));

hits_data = ones(size(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Mat_a_obs_v2.0.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

ralser_lehrach_2008.orfs = tested_orfs;
ralser_lehrach_2008.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
ralser_lehrach_2008.data(ind2,:) = hits_data(ind1,:);

ralser_lehrach_2008.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'ralser_lehrach_2008.mat'],'ralser_lehrach_2008');
return;

% Save data into database
dt = ralser_lehrach_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

