%% Teixeira~Sa-Correia, 2009
% DATA = teixeira_sa_correia_2009
function FILENAMES = code()
FILENAMES = {};
teixeira_sa_correia_2009.pmid = 19633105;

phenotypes = {'growth'};
treatments = {'EtOH'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','/Users/brianna/Documents/Datasets/Phenotypes/2009_Teixeira~Sa-Correia/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','/Users/brianna/Documents/Datasets/Phenotypes/2009_Teixeira~Sa-Correia/TableS1_suplementary_material.xlsx');
hits_genenames = data.raw(7:end,1);
inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];

hits_orfs = genename2orf(hits_genenames);
hits_orfs(strcmp('opi8', hits_orfs))= {'YKR035C'};
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('ppa1', hits_orfs)) = [];  % ambiguous genename
hits_orfs(strcmp('rlr1', hits_orfs)) = {'YNL139C'};
hits_orfs(strcmp('soy1', hits_orfs)) = {'YBR194W'};
hits_orfs(strcmp('tfp3', hits_orfs)) = {'YPL234C'};
hits_orfs(strcmp('vps66', hits_orfs)) = {'YPR139C'};

hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];


[missing, ix] = setdiff(hits_orfs, tested_orfs);
hits_orfs(ix) = [];      % 24 ORFs eliminated from the hit list (the list contains hits from both the hap and the het collection, but the list of tested strains is only available for the hap collection)

teixeira_sa_correia_2009.orfs = tested_orfs;
teixeira_sa_correia_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
teixeira_sa_correia_2009.data(ind2) = -1;

teixeira_sa_correia_2009.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'teixeira_sa_correia_2009.mat'],'teixeira_sa_correia_2009');
return;

% Save data into database
dt = teixeira_sa_correia_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
database_ix = database_ix(2);
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

