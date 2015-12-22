%% Burston~Conibear, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
burston_conibear_2009.pmid = 19506040;

phenotypes = {'endocytosis (MatA)';'endocytosis (MatAlpha)'};
treatments = {''};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/JCB_200811116_TS1.xlsx', 'TableS1');

hits_orfs = data.raw(6:end,2);
hits_data_a = data.raw(6:end,4);
hits_data_alpha = data.raw(6:end,5);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

burston_conibear_2009.orfs = hits_orfs;
burston_conibear_2009.data = cell2mat([hits_data_a hits_data_alpha]);

burston_conibear_2009.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'burston_conibear_2009.mat'],'burston_conibear_2009');
return;

% Save data into database
dt = burston_conibear_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

