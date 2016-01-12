%% Merz~Westermann, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
merz_westermann_2009.pmid = 19751518;

% Part 1
phenotypes = {'growth'};
treatments = {'Gly, 3%'};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/gb-2009-10-9-r95-s1.xlsx');
hits_orfs = data.raw(3:end,1);
hits_orfs = unique(strtrim(upper(hits_orfs)));

hits_data = -ones(size(hits_orfs));

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/pet-Screen.xlsx', 'mat_alpha_obs');
tested_orfs = tested.raw(2:end,2);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

merz_westermann_2009.orfs = tested_orfs;
merz_westermann_2009.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
merz_westermann_2009.data(ind2,:) = hits_data(ind1,:);

merz_westermann_2009.ph = strcat(phenotypes, '; ', treatments);

save('./merz_westermann_2009.mat','merz_westermann_2009');
return;

% Save data into database
dt = merz_westermann_2009;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

