%% Copic~Miller, 2009
% DATA = copic_miller_2009
function FILENAMES = code()
FILENAMES = {};
copic_miller_2009.pmid = 19433630;

phenotypes = {'Kar2 secretion'};
treatments = {''};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/mat_alpha_obs_v1.0.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,1);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/TableS1.xlsx', 'TABLE S1');
hits_orfs = data.raw(5:end,1);
hits_notes = data.raw(5:end,3);
hits_scores = data.raw(5:end,6);

inds = strcmp('mat-a only', hits_notes);
hits_orfs(inds) = [];
hits_scores(inds) = [];

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds) = [];
hits_orfs = upper(hits_orfs);

hits_scores = cell2mat(hits_scores);

copic_miller_2009.orfs = tested_orfs;
copic_miller_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
copic_miller_2009.data(ind2,1) = hits_scores(ind1);

copic_miller_2009.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'copic_miller_2009.mat'],'copic_miller_2009');
return;

% Save data into database
dt = copic_miller_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

