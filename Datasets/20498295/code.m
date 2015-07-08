%% Dakshinamurthy~Garfinkel, 2010
% DATA = dakshinamurthy_garfinkel_2010
function FILENAMES = code()
FILENAMES = {};
dakshinamurthy_garfinkel_2010.pmid = 20498295;

phenotypes = {'decreased Ty1 transposon mobility'};
treatments = {''};

load Datasets/Phenotypes_standard/18202368/nyswaner_garfinkel_2008;
dakshinamurthy_garfinkel_2010.orfs = nyswaner_garfinkel_2008.orfs;      % Same screen, positive & negative phenotypes analyzed separately
dakshinamurthy_garfinkel_2010.data = zeros(length(dakshinamurthy_garfinkel_2010.orfs),1);

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/TableS1-2.xlsx', 'Sheet2');
hits_orfs = data.raw(:,2);
hits_data = data.raw(:,3);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_data = cellfun(@length, hits_data);

missing = setdiff(hits_orfs, dakshinamurthy_garfinkel_2010.orfs);

% Adding 4 ORFs to the list of tested
dakshinamurthy_garfinkel_2010.orfs = [dakshinamurthy_garfinkel_2010.orfs; missing];

[~,ind1,ind2] = intersect(dakshinamurthy_garfinkel_2010.orfs, hits_orfs);
dakshinamurthy_garfinkel_2010.data(ind1,1) = hits_data(ind2);

dakshinamurthy_garfinkel_2010.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'dakshinamurthy_garfinkel_2010.mat'],'dakshinamurthy_garfinkel_2010');
return;

% Save data into database
dt = dakshinamurthy_garfinkel_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

