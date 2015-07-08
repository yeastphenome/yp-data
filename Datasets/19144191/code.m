%% Kemmer~Roberge, 2009
% DATA = kemmer_roberge_2009
function FILENAMES = code()
FILENAMES = {};

kemmer_roberge_2009.pmid = 19144191;

phenotypes = {'growth (colony size)'};
treatments = {'dhMotC [60 uM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/haploid set.xlsx', 'haploid set');
tested_orfs = tested.raw(6:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','raw_data/hits_genenames.txt', '%s');

hits_genenames = cellfun(@strtrim, hits_genenames,'UniformOutput',0);
hits_orfs = genename2orf(hits_genenames,'noannot');
hits_orfs = unique(upper(hits_orfs));

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

kemmer_roberge_2009.orfs = tested_orfs;
kemmer_roberge_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
kemmer_roberge_2009.data(ind2,1) = hits_scores(ind1);

kemmer_roberge_2009.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'kemmer_roberge_2009.mat'],'kemmer_roberge_2009');
return;

% Save data into database
dt = kemmer_roberge_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

