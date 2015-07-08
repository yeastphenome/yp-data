%% Begley~Samson, 2004
% DATA = begley_samson_2004
function FILENAMES = code()
FILENAMES = {};
begley_samson_2004.pmid = 15469827;

phenotypes = {'growth (spot assay)'};
treatments = {'MMS';'t-BuOOH';'4NQO';'UV score'};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/Begley2003.xlsx', 'Sheet1');

hits_orfs = data.raw(2:end,1);
hits_data = data.raw(2:end,7:10);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_orfs = upper(cellfun(@strtrim, hits_orfs,'UniformOutput',0));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = -cell2mat(hits_data);   % The sensitivity scores are such that the higher the number, the more sensitive the mutant.

[t,t2] = grpstats(hits_data,hits_orfs,{'gname','mean'});
begley_samson_2004.orfs = t;
begley_samson_2004.data = t2;

begley_samson_2004.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'begley_samson_2004.mat'],'begley_samson_2004');
return;

% Save data into database
dt = begley_samson_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

