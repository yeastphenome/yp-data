%% Teixeira~Sa-Correia, 2010
% DATA = teixeira_sa_correia_2010
function FILENAMES = code()
FILENAMES = {};
teixeira_sa_correia_2010.pmid = 20210661;

phenotypes = {'growth'};
treatments = {'Glu, 30%'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','/Users/brianna/Documents/Datasets/Phenotypes/2010_Teixeira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, hits_genenames] = dataread('textread','/Users/brianna/Documents/Datasets/Phenotypes/2010_Teixeira~Sa-Correia/hits_genenames.txt', '%s');
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_genenames = strtrim(hits_genenames);

hits_orfs = genename2orf(hits_genenames);
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('ppa1', hits_orfs)) = [];  % ambiguous genename
hits_orfs(strcmp('soy1', hits_orfs)) = {'YBR194W'};

hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];

hits_orfs = unique(hits_orfs);


[missing, ix] = setdiff(hits_orfs, tested_orfs);
hits_orfs(strcmp('YML009W-B',hits_orfs)) = {'YML010W-A'};

teixeira_sa_correia_2010.orfs = tested_orfs;
teixeira_sa_correia_2010.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2010.data(ind2) = -1;

teixeira_sa_correia_2010.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'teixeira_sa_correia_2010.mat'],'teixeira_sa_correia_2010');
return;

% Save data into database
dt = teixeira_sa_correia_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

