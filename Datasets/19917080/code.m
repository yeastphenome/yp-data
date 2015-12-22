%% Arita~Costa, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
arita_costa_2009.pmid = 19917080;

phenotypes = {'growth [spot assay]'};
treatments = {'NiSO4 [0.75-1.25 mM]'};

% Load resistant
[FILENAMES{end+1}, hits.raw] = dataread('xlsread','./raw_data/1471-2164-10-524-S1.XLS', 'Sheet1');
hits_sensitive_orfs = hits.raw(7:end,1);
hits_resistant_orfs = hits.raw(7:end,2);

inds = cellfun(@isnumeric, hits_resistant_orfs);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = cellfun(@strtrim, hits_resistant_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_resistant_orfs,1);
hits_resistant_orfs(inds) = [];
hits_resistant_orfs = unique(upper(hits_resistant_orfs));

hits_resistant_scores = ones(length(hits_resistant_orfs),1);

inds = cellfun(@isnumeric, hits_sensitive_orfs);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = cellfun(@strtrim, hits_sensitive_orfs,'UniformOutput',0);
inds = ~strncmp('Y', hits_sensitive_orfs,1);
hits_sensitive_orfs(inds) = [];
hits_sensitive_orfs = unique(upper(hits_sensitive_orfs));

hits_sensitive_scores = -ones(length(hits_sensitive_orfs),1);



% Check overlap between resistant and sensitive
length(intersect(hits_resistant_orfs, hits_sensitive_orfs))

arita_costa_2009.orfs = unique([hits_resistant_orfs; hits_sensitive_orfs]);
arita_costa_2009.data = zeros(length(arita_costa_2009.orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_resistant_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2,1) = hits_resistant_scores(ind1);
[~,ind1,ind2] = intersect(hits_sensitive_orfs, arita_costa_2009.orfs);
arita_costa_2009.data(ind2,1) = hits_sensitive_scores(ind1);


arita_costa_2009.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'arita_costa_2009.mat'],'arita_costa_2009');
return;

% Save data into database
dt = arita_costa_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

