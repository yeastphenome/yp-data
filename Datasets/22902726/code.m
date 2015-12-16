%% Schmidt~Boyer, 2012
function FILENAMES = code()
FILENAMES = {};

schmidt_boyer_2012.pmid = 22902726;

phenotypes = {'growth (OD)';'growth (colony size)';'growth (MIC)'};
treatments = {'boric acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/BA sensitivity comprehensive data file.xlsx', 'Strain list');
tested_orfs = tested.raw(2:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/BA sensitivity comprehensive data file.xlsx', 'Sheet1');

% Dataset1: Tested = all; hits = liquid
hits_orfs = data.raw(3:end,1);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));
hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

schmidt_boyer_2012.orfs = tested_orfs;
schmidt_boyer_2012.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,1) = hits_scores(ind1);

% Dataset2: Tested = dataset1, hits = solid
hits_orfs = data.raw(3:end,3);
inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_orfs = unique(upper(hits_orfs));
hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(schmidt_boyer_2012.data(:,1)==0);
schmidt_boyer_2012.data(inds,2) = NaN;
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,2) = hits_scores(ind1);


% Dataset3: Tested = dataset2, data = MIC
hits_orfs = data.raw(3:end,5);
hits_scores = data.raw(3:end,7);

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

hits_orfs = unique(upper(hits_orfs));
hits_scores = -cell2mat(hits_scores); %MIC values transformed to -MIC such that lower values correspond to the most sensitive strains (high original MIC)

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(schmidt_boyer_2012.data(:,2)==0 | isnan(schmidt_boyer_2012.data(:,2)));
schmidt_boyer_2012.data(inds,3) = NaN;
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
schmidt_boyer_2012.data(ind2,3) = hits_scores(ind1);

schmidt_boyer_2012.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'schmidt_boyer_2012.mat'],'schmidt_boyer_2012');
return;


% Save data into database
dt = schmidt_boyer_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
adj_ix = [1 3 2];
[~, adj_ix] = sort(adj_ix);

datasets.names(database_ix(adj_ix),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix(adj_ix)));

end

