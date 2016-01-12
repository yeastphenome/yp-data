%% Thevissen~Francois, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

thevissen_francois_2007.pmid = 17553796;

phenotypes = {'growth [MIC]'};
treatments = {'miconazole [0.025-12.5 ug/ml]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/Euroscarf library.xlsx', 'Tabelle1');
tested_orfs = tested.raw(2:end,2);


inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = ~strncmp('Y', tested_orfs,1);
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
fid = fopen('raw_data/hits_orfs_scores.txt');
hits = textscan(fid,'%s %d');
hits_orfs = upper(hits{1});
hits_scores = -hits{2};
fclose(fid);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustment
tested_orfs = [tested_orfs; missing];   % 2 ORFs added

thevissen_francois_2007.orfs = tested_orfs;
thevissen_francois_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
thevissen_francois_2007.data(ind2,1) = hits_scores(ind1);

thevissen_francois_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('./thevissen_francois_2007.mat','thevissen_francois_2007');
return;

% Save data into database
dt = thevissen_francois_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

