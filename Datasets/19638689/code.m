%% Auesukaree~Harashima, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
auesukaree_harashima_2009.pmid = 19638689;

phenotypes = {'growth'};
treatments = {'EtOH';'MeOH';'propanol';'NaCl';'H2O2';'37C'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Mat alpha_KOset list.xlsx');
tested_orfs = tested.raw(4:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = upper(regexprep(tested_orfs, '\W',''));
inds = find(~isorf(tested_orfs));
for i = 1 : length(inds)
    tested_orfs{inds(i)} = [tested_orfs{inds(i)}(1:end-1) '-' tested_orfs{inds(i)}(end)];
end

% Load data
[FILENAMES{end+1}, data_hits{1}] = dataread('textread','./raw_data/ethanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{2}] = dataread('textread','./raw_data/methanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{3}] = dataread('textread','./raw_data/propanol_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{4}] = dataread('textread','./raw_data/nacl_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{5}] = dataread('textread','./raw_data/h2o2_sensitivity_hits.txt', '%s');
[FILENAMES{end+1}, data_hits{6}] = dataread('textread','./raw_data/heat_sensitivity_hits.txt', '%s');

data_hits_orfs = cell(size(data_hits));

data_hits_orfs{1} = unique(translate(data_hits{1}));
[missing, ix] = setdiff(data_hits_orfs{1}, tested_orfs);
data_hits_orfs{1}(strcmp('YHR039C-A', data_hits_orfs{1})) = {'YHR039C-B'};

data_hits_orfs{2} = unique(translate(data_hits{2}));
[missing, ix] = setdiff(data_hits_orfs{2}, tested_orfs);

data_hits_orfs{3} = unique(translate(data_hits{3}));
[missing, ix] = setdiff(data_hits_orfs{3}, tested_orfs);

data_hits_orfs{4} = unique(translate(data_hits{4}));
[missing, ix] = setdiff(data_hits_orfs{4}, tested_orfs);

data_hits_orfs{5} = unique(translate(data_hits{5}));
[missing, ix] = setdiff(data_hits_orfs{5}, tested_orfs);


data_hits_orfs{6} = unique(translate(data_hits{6}));
[missing, ix] = setdiff(data_hits_orfs{6}, tested_orfs);
tested_orfs = [tested_orfs; missing];

auesukaree_harashima_2009.orfs = tested_orfs;
auesukaree_harashima_2009.data = zeros(length(tested_orfs),length(treatments));

for i = 1 : length(treatments)
    [~,ind1,ind2] = intersect(data_hits_orfs{i}, tested_orfs);
    auesukaree_harashima_2009.data(ind2,i) = -1;
end

auesukaree_harashima_2009.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'auesukaree_harashima_2009.mat'],'auesukaree_harashima_2009');
return;

% Save data into database
dt = auesukaree_harashima_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
ph_ix = ph_ix([2 3 4 5 6 1]);
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

