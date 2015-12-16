%% Bishop~Avery, 2007
function FILENAMES = code()
FILENAMES = {};

bishop_avery_2007.pmid = 17176259;

phenotypes = {'growth (colony size)'};
treatments = {'Ni(NO3)2 [2.5-4 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Mata_DeletionArray+slow_growers.xlsx', '96');
tested_orfs = tested.raw(3:end,2);
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Mata_DeletionArray+slow_growers.xlsx', 'slow growers');
tested_orfs = [tested_orfs; tested.raw(3:end,2)];
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, hits_genenames_resistant] = dataread('textread','./raw_data/hits_genenames_resistant.txt', '%s');

hits_genenames_resistant = upper(hits_genenames_resistant);
hits_orfs_resistant = genename2orf(hits_genenames_resistant,'noannot');
hits_orfs_resistant = cellfun(@strtrim, hits_orfs_resistant,'UniformOutput',0);
hits_orfs_resistant = unique(upper(hits_orfs_resistant));

hits_scores_resistant = ones(length(hits_orfs_resistant),1);

inds = find(~strncmp('Y', hits_orfs_resistant,1));
hits_orfs_resistant(inds) = [];
hits_scores_resistant(inds) = [];


[FILENAMES{end+1}, hits_genenames_sensitive] = dataread('textread','./raw_data/hits_genenames_sensitive.txt', '%s');

hits_genenames_sensitive = upper(hits_genenames_sensitive);
hits_orfs_sensitive = genename2orf(hits_genenames_sensitive,'noannot');

% Adjustments
hits_orfs_sensitive(strcmpi('FMP18', hits_orfs_sensitive)) = {'YKR065C'};

hits_orfs_sensitive = cellfun(@strtrim, hits_orfs_sensitive,'UniformOutput',0);
hits_orfs_sensitive = unique(upper(hits_orfs_sensitive));

hits_scores_sensitive = -ones(length(hits_orfs_sensitive),1);

inds = find(~strncmp('Y', hits_orfs_sensitive,1));
hits_orfs_sensitive(inds) = [];
hits_scores_sensitive(inds) = [];

overlapping = intersect(hits_orfs_resistant, hits_orfs_sensitive);

hits_orfs = [hits_orfs_resistant; hits_orfs_sensitive];
hits_scores = [hits_scores_resistant; hits_scores_sensitive];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 5 orfs to be added

bishop_avery_2007.orfs = tested_orfs;
bishop_avery_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
bishop_avery_2007.data(ind2,1) = hits_scores(ind1);

bishop_avery_2007.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'bishop_avery_2007.mat'],'bishop_avery_2007');
return;

% Save data into database
dt = bishop_avery_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

