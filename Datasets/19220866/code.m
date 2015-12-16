%% Mira~Sa-Correia, 2009
function FILENAMES = code()
FILENAMES = {};
mira_sa_correia_2009.pmid = 19220866;

phenotypes = {'growth'};
treatments = {'propionic acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','./raw_data/hits_genenames.txt', '%s');
[FILENAMES{end+1}, hits_genenames_moderate] = dataread('textread','./raw_data/hits_genenames_moderate.txt', '%s');

hits_genenames = strtrim(lower(hits_genenames));
hits_genenames_moderate = strtrim(lower(hits_genenames_moderate));

hits_genenames_strong = setdiff(hits_genenames, hits_genenames_moderate);


hits_orfs_moderate = genename2orf(hits_genenames_moderate);
hits_orfs_moderate(strcmp('set7', hits_orfs_moderate)) = {'YDR257C'};
hits_orfs_moderate = unique(upper(hits_orfs_moderate));

hits_orfs_strong = genename2orf(hits_genenames_strong);
hits_orfs_strong(strcmp('mrp16', hits_orfs_strong)) = [];    % typo? genename doesn't exist
hits_orfs_strong(strcmp('rhr2', hits_orfs_strong)) = {'YIL053W'};
hits_orfs_strong(strcmp('tfp1', hits_orfs_strong)) = {'YDL185W'};
hits_orfs_strong(strcmp('tfp3', hits_orfs_strong)) = {'YPL234C'};
hits_orfs_strong = unique(upper(hits_orfs_strong));

inds = find(~strncmp('Y', hits_orfs_moderate,1));
hits_orfs_moderate(inds) = [];

inds = find(~strncmp('Y', hits_orfs_strong,1));
hits_orfs_strong(inds) = [];

[missing, ix] = setdiff(hits_orfs_strong, tested_orfs);
hits_orfs_strong(ix) = [];      % 3 ORFs eliminated from the hit list

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

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'mira_sa_correia_2009.mat'],'mira_sa_correia_2009');
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

