%% Postma~Ralser, 2009
% DATA = postma_ralser_2009
function FILENAMES = code()
FILENAMES = {};

postma_ralser_2009.pmid = 20157578;

phenotypes = {'growth [spot assay]'};
treatments = {'hibernation'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/Mat_a_obs_v2.0.xlsx', 'Mat_a_obs_v2.0.txt');
tested_orfs = tested.raw(2:end,1);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('raw_data/hits_orfs.txt');
hits_orfs = textscan(fid,'%s');
fclose(fid);

hits_orfs = unique(hits_orfs{1});
hits_scores = zeros(length(hits_orfs),1)+1;

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 1 ORF added

postma_ralser_2009.orfs = tested_orfs;
postma_ralser_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
postma_ralser_2009.data(ind2,1) = hits_scores(ind1);

postma_ralser_2009.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'postma_ralser_2009.mat'],'postma_ralser_2009');
return;

% Save data into database
dt = postma_ralser_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

