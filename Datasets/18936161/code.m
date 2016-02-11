%% Mir~Cashikar, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

mir_cashikar_2009.pmid = 18936161;

phenotypes = {'growth [MIC]'};
treatments = {'heat stress (temperature [50ºC], duration [30 min])'};

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/Mat_a.xlsx', 'mat_a_041902');
tested_orfs = tested.raw(3:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('./raw_data/hits_orfs_scores.txt');
hits = textscan(fid,'%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};
fclose(fid);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

mir_cashikar_2009.orfs = tested_orfs;
mir_cashikar_2009.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mir_cashikar_2009.data(ind2,1) = hits_scores(ind1);

mir_cashikar_2009.ph = [strcat(phenotypes, '; ', treatments)];

save('./mir_cashikar_2009.mat','mir_cashikar_2009');
return;

% Save data into database
dt = mir_cashikar_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

