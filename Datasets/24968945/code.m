%% Gaupel~Tenniswood, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gaupel_tenniswood_2004.pmid = 24968945;

phenotypes = {'growth'};
treatments = {'CG-1521'};

% Load data (1)
[FILENAMES{end+1}, data_sens.raw] = readdata('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top sensitive strains');

hits_orfs1 = data_sens.raw(2:end,1);
hits_data1 = -cell2mat(data_sens.raw(2:end,3));

hits_orfs1 = strtrim(upper(hits_orfs1));

inds = find(~strncmp('Y', hits_orfs1,1));
hits_orfs1(inds) = [];
hits_data1(inds,:) = [];

[t,t2] = grpstats(hits_data1,hits_orfs1,{'gname','mean'});
hits_orfs1 = t;
hits_data1 = t2;

% Load data (2)
[FILENAMES{end+1}, data_res.raw] = readdata('xlsread','./raw_data/1471-2164-15-528-s1.xlsx', 'Top resistant strains');

hits_orfs2 = data_res.raw(2:end,1);
hits_data2 = -cell2mat(data_res.raw(2:end,3));

hits_orfs2 = strtrim(upper(hits_orfs2));

inds = find(~strncmp('Y', hits_orfs2,1));
hits_orfs2(inds) = [];
hits_data2(inds,:) = [];

[t,t2] = grpstats(hits_data2,hits_orfs2,{'gname','mean'});
hits_orfs2 = t;
hits_data2 = t2;

% Eliminate overlap between sensitive and resistant strains?
[~,ind1,ind2] = intersect(hits_orfs1, hits_orfs2);
hits_orfs1(ind1) = []; hits_data1(ind1,:) = [];
hits_orfs2(ind2) = []; hits_data2(ind2,:) = [];

% Combine part1 and part2
hits_orfs = [hits_orfs1; hits_orfs2];
hits_data = [hits_data1; hits_data2];

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/CompleteDeletionLibrary.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

gaupel_tenniswood_2004.orfs = tested_orfs;
gaupel_tenniswood_2004.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
gaupel_tenniswood_2004.data(ind2,:) = hits_data(ind1,:);

gaupel_tenniswood_2004.ph = strcat(phenotypes, '; ', treatments);

save('./gaupel_tenniswood_2004.mat','gaupel_tenniswood_2004');
return;

% Save data into database
dt = gaupel_tenniswood_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

