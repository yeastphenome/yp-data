%% Mira~Sa-Correia, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mira_sa_correia_2009.pmid = 19220866;

phenotypes = {'growth'};
treatments = {'propionic acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
[FILENAMES{end+1}, hits_genenames] = readdata('textread','./raw_data/hits_genenames.txt', '%s');
[FILENAMES{end+1}, hits_genenames_moderate] = readdata('textread','./raw_data/hits_genenames_moderate.txt', '%s');

hits_genenames = cleanGenename(hits_genenames);
hits_genenames_moderate = cleanGenename(hits_genenames_moderate);

hits_genenames_strong = setdiff(hits_genenames, hits_genenames_moderate);


hits_orfs_moderate = translate(hits_genenames_moderate);
hits_orfs_moderate = unique(hits_orfs_moderate);

[hits_orfs_strong, translated] = translate(hits_genenames_strong);
hits_orfs_strong(~translated) = [];
hits_orfs_strong = unique(hits_orfs_strong);

[missing, ix] = setdiff(hits_orfs_strong, tested_orfs);
hits_orfs_strong(ix) = [];      % 1 ORF eliminated from the hit list

[missing, ix] = setdiff(hits_orfs_moderate, tested_orfs);

hits_data_strong = zeros(length(hits_orfs_strong),1)-2;
hits_data_moderate = zeros(length(hits_orfs_moderate),1)-1;


mira_sa_correia_2009.orfs = tested_orfs;
mira_sa_correia_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs_strong, tested_orfs);
mira_sa_correia_2009.data(ind2) = hits_data_strong(ind1);

[~,ind1,ind2] = intersect(hits_orfs_moderate, tested_orfs);
mira_sa_correia_2009.data(ind2) = hits_data_moderate(ind1);


mira_sa_correia_2009.ph = strcat(phenotypes, '; ', treatments);

save('./mira_sa_correia_2009.mat','mira_sa_correia_2009');
return;

% Save data into database
dt = mira_sa_correia_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

