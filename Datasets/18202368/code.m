%% Nyswaner~Garfinkel,2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
nyswaner_garfinkel_2008.pmid = 18202368;

phenotypes = {'increased Ty1 transposon mobility'};
treatments = {''};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Matalphakos counted.xlsx', 'Sheet2');
tested_orfs = tested.raw(6:end,2);

inds = find(cellfun(@isempty, tested_orfs));
tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/genetics.107.082602-9.xlsx', 'Sheet1');
hits_orfs = data.raw(:,1);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];

nyswaner_garfinkel_2008.orfs = tested_orfs;
nyswaner_garfinkel_2008.data = zeros(length(tested_orfs),1);

missing = setdiff(hits_orfs, tested_orfs);

% Adding 3 ORFs to the list of tested
tested_orfs = [tested_orfs; missing];

[~,ind1,ind2] = intersect(tested_orfs, hits_orfs);
nyswaner_garfinkel_2008.data(ind1,1) = 1;

nyswaner_garfinkel_2008.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'nyswaner_garfinkel_2008.mat'],'nyswaner_garfinkel_2008');
return;

% Save data into database
dt = nyswaner_garfinkel_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

