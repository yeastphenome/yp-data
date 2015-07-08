%% Brett~Rao, 2011
% DATA = brett_rao_2011
function FILENAMES = code()
FILENAMES = {};

brett_rao_2011.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0017619.s003'};
brett_rao_2011.downloaddate = {'2013-07-18'};
brett_rao_2011.pmid = 21423800;

phenotypes = {'Growth, OD600';'Vacuolar pH (pHv)'};
treatments = {'external pH 2.7';'external pH 4.0';'external pH 7.0'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/journal.pone.0017619.s003.xlsx', 'Unsorted Data');

% Get indices of the data columns
ind_data = [4:6 7:9];


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

inds = find(~strncmp('Y', data.raw(:,2),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,2) = cellfun(@strtrim, data.raw(:,2),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,2));
data2.data = cell2mat(t);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
brett_rao_2011.orfs = t;
brett_rao_2011.data = t2;
brett_rao_2011.ph = [strcat(phenotypes{1}, '; ', treatments); strcat(phenotypes{2},'; ', treatments)];


% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(brett_rao_2011.orfs, essential_genes);
brett_rao_2011.orfs(ind1) = [];
brett_rao_2011.data(ind1,:) = [];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'brett_rao_2011.mat'],'brett_rao_2011');
return;

% Save data into database
dt = brett_rao_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

