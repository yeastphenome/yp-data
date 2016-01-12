%% Chamilos~Kontoyiannis, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chamilos_kontoyiannis_2008.pmid = 18212113;

phenotypes = {'growth'};
treatments = {'gliotoxin'};

% Load tested
[FILENAMES{end+1}, tested_orfs] = dataread('textread','./raw_data/tested_orfs.txt', '%s');

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
fid = fopen('./raw_data/data_genenames.txt','r');
C = textscan(fid, '%s %f','delimiter','\t');
fclose(fid);

genenames = C{1};
raw_data = C{2};

% Normalize to WT
inds = find(strcmp('WT', genenames));
raw_data = raw_data./raw_data(inds)-1;

genenames(inds) = [];
raw_data(inds) = [];

orfs = translate(genenames);

missing_orfs = setdiff(orfs, tested_orfs);
tested_orfs = [tested_orfs; missing_orfs];  % Adding 2 missing strains to the list of tested strains.

chamilos_kontoyiannis_2008.orfs = tested_orfs;
chamilos_kontoyiannis_2008.data = zeros(length(tested_orfs),1);
[~,ind1,ind2] = intersect(tested_orfs, orfs);
chamilos_kontoyiannis_2008.data(ind1) = raw_data(ind2);
chamilos_kontoyiannis_2008.ph = strcat(phenotypes, '; ', treatments);

save('./chamilos_kontoyiannis_2008.mat','chamilos_kontoyiannis_2008');
return;

% Save data into database
dt = chamilos_kontoyiannis_2008;

datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

