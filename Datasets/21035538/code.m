%% Uluisik~Koc, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
uluisik_koc_2011.pmid = 21035538;

phenotypes = {'growth'};
treatments = {'BA, 20 mM';'BA, 30 mM';'BA, 40 mM';'BA, 50 mM'; 'BA, 60 mM'; 'BA, 70 mM';'BA, 100 mM';'BA, 125 mM';'BA, 150 mM'};

% Load tested
[FILENAMES{end+1}, tested_orfs] = readdata('textread','./raw_data/tested_strains.txt', '%s');

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(strtrim(tested_orfs)));

uluisik_koc_2011.orfs = tested_orfs;
uluisik_koc_2011.data = zeros(length(tested_orfs),length(treatments));
uluisik_koc_2011.ph = strcat(phenotypes, '; ', treatments);

% Load data 1
[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/mmc2.xlsx');
hits_orfs = data.raw(4:end,1);
hits_data = data.raw(4:end,2:5);
hits_treatments = {'BA, 100 mM','BA, 125 mM','BA, 150 mM'};

hits_data(strcmp('+',hits_data)) = {1};
hits_data(strcmp('-',hits_data)) = {0};

hits_data = cell2mat(hits_data);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = hits_data - repmat(hits_data(1,:),size(hits_data,1),1);    % Subtract the 1st row (WT)
hits_data = hits_data - repmat(hits_data(:,1),1,size(hits_data,2));     % Subtract the 1st col (0 mM);
hits_data(1,:) = [];
hits_data(:,1) = [];
hits_orfs(1) = [];

hits_orfs = strtrim(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];



[missing, ix] = setdiff(hits_orfs, tested_orfs);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
[~,ind3,ind4] = intersect(hits_treatments, treatments);
uluisik_koc_2011.data(ind2,ind4) = hits_data(ind1,ind3);

% Load data 2
[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/mmc3.xlsx');
hits_orfs = data.raw(4:end,1);
hits_data = data.raw(4:end,2:8);
hits_treatments = {'BA, 20 mM';'BA, 30 mM';'BA, 40 mM';'BA, 50 mM';'BA, 60 mM';'BA, 70 mM'};

hits_data(strcmp('+',hits_data)) = {1};
hits_data(strcmp('-',hits_data)) = {0};

hits_data = cell2mat(hits_data);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data = hits_data - repmat(hits_data(1,:),size(hits_data,1),1);    % Subtract the 1st row (WT)
hits_data = hits_data - repmat(hits_data(:,1),1,size(hits_data,2));     % Subtract the 1st col (0 mM);
hits_data(1,:) = [];
hits_data(:,1) = [];
hits_orfs(1) = [];

hits_orfs = strtrim(upper(hits_orfs));
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 3 ORFs missing (eliminated)
hits_orfs(ix) = [];
hits_data(ix,:) = [];

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
[~,ind3,ind4] = intersect(hits_treatments, treatments);
uluisik_koc_2011.data(ind2,ind4) = hits_data(ind1,ind3);


save('./uluisik_koc_2011.mat','uluisik_koc_2011');
return;

% Save data into database
dt = uluisik_koc_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

