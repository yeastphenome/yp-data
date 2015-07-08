%% Orij~Smits, 2012
% DATA = orij_smits_2012
function FILENAMES = code()
FILENAMES = {};

orij_smits_2012.pmid = 23021432;

phenotypes = {'cytosolic pH'};
treatments = {'standard'};

% Load hits
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/pH Screen raw.xlsx', 'initial screens');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3:4);

inds = cellfun(@isnumeric, hits_orfs);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_orfs,1);
hits_orfs(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = upper(hits_orfs);
hits_scores = cell2mat(hits_scores);
hits_scores = nanmean(hits_scores, 2);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

orij_smits_2012.orfs = hits_orfs;
orij_smits_2012.data = hits_scores;

orij_smits_2012.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'orij_smits_2012.mat'],'orij_smits_2012');
return;

% Save data into database
dt = orij_smits_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

