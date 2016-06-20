%% Tigano~Ottonello, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
tigano_ottonello_2015.pmid = 26240381;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(tigano_ottonello_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load hit strains and tested strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/nar-03484-z-2014-File007.xlsx', 'Suppl_TabS1');
[FILENAMES{end+1}, dataAll] = read_data('xlsread', './raw_data/Tigano et al NAR 2015_List of tested yeast mutant strains.xlsx', 'Foglio1');

% Get list of ORFs from both sets of ORFs
all_orfs = dataAll(2:end, 1);
hit_orfs = data(4:end,1);
hit_data = data(4:end, 8:9);

% Get rid of empty rows
inds = find(cellfun(@isnumeric, all_orfs));
all_orfs(inds) = [];

% Clean up ORFs
hit_orfs = clean_orf(hit_orfs);
all_orfs = clean_orf(all_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds)); 
inds = find(~is_orf(all_orfs));
disp(all_orfs(inds));

% Make sure all are unique
all_orfs = unique(all_orfs);

% Set Numeric Values to hits data
hit_data(strcmp('HS', hit_data)) = {-3};
hit_data(strcmp('MS', hit_data)) = {-2};
hit_data(strcmp('LS', hit_data)) = {-1};

% Transform cell array into Double Array
hit_data = cell2mat(hit_data);

% Average any repeated value
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% Check to see all hit_orfs are in all_orfs
missing_orfs = setdiff(hit_orfs, all_orfs);     % 1 ORF missing
all_orfs = [all_orfs; missing_orfs];

% Make a zero matrix for all data points
all_data = zeros(length(all_orfs), 2);

% Find indices for hit_orfs in all_orfs 
[~, ind1, ind2] = intersect(all_orfs, hit_orfs);
all_data(ind1,:) = hit_data(ind2,:);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [695; 694];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

tigano_ottonello_2015.orfs = all_orfs;
tigano_ottonello_2015.ph = hit_data_names;
tigano_ottonello_2015.data = all_data;
tigano_ottonello_2015.dataset_ids = hit_data_ids;

%% Save

save('./tigano_ottonello_2015.mat','tigano_ottonello_2015');

fid = fopen('./tigano_ottonello_2015.txt','w');
write_matrix_file(fid, tigano_ottonello_2015.orfs, tigano_ottonello_2015.ph, tigano_ottonello_2015.data);
fclose(fid);

end
