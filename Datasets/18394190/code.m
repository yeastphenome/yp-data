%% Ruotolo~Ottonello, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ruotolo_ottonello_2008.pmid = 18394190;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ruotolo_ottonello_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/gb-2008-9-4-r67-S2.xlsx', 'Additional file 2');

% Get the list of ORFs and the correponding data 
hit_strains = data(4:end,1);

% Get the data itself
hit_data = cell2mat(data(4:end,[2 4]));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1317; 1318];

%% Load tested strains

[FILENAMES{end+1}, dataAll] = read_data('xlsread', './raw_data/List of tested yeast mutant strains.xlsx', 'Foglio1');

all_orfs = dataAll(2:end, 1);

% Get rid of empty rows
inds = find(cellfun(@isnumeric, all_orfs));
all_orfs(inds) = [];

% Clean up ORFs
all_orfs = clean_orf(all_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(all_orfs));
disp(all_orfs(inds));

% Make sure all are unique
all_orfs = unique(all_orfs);

% Check to see all hit_orfs are in all_orfs
missing_orfs = setdiff(hit_strains, all_orfs);     % 5 ORFs missing
all_orfs = [all_orfs; missing_orfs];

% Make a zero matrix for all data points
all_data = zeros(length(all_orfs), 2);

% Find indices for hit_orfs in all_orfs 
[~, ind1, ind2] = intersect(all_orfs, hit_strains);
all_data(ind1,:) = hit_data(ind2,:);

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ruotolo_ottonello_2008.orfs = all_orfs;
ruotolo_ottonello_2008.ph = hit_data_names;
ruotolo_ottonello_2008.data = all_data;
ruotolo_ottonello_2008.dataset_ids = hit_data_ids;

%% Save

save('./ruotolo_ottonello_2008.mat','ruotolo_ottonello_2008');

%% Print out

fid = fopen('./ruotolo_ottonello_2008.txt','w');
write_matrix_file(fid, ruotolo_ottonello_2008.orfs, ruotolo_ottonello_2008.ph, ruotolo_ottonello_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ruotolo_ottonello_2008)
end

end

