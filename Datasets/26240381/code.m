%% Tigano~Ottonello, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
tigano_ottonello_2015.pmid = 26240381;

phenotypes = {'growth'};
treatments = {'ethanol, 30C', 'ethanol, 37C'};

%% Hit Strains

% Load hit strains and tested strains
[FILENAMES{end+1}, data] = readdata('xlsread', './raw_data/nar-03484-z-2014-File007.xlsx', 'Suppl_TabS1');
[FILENAMES{end+1}, dataAll] = readdata('xlsread', './raw_data/Tigano et al NAR 2015_List of tested yeast mutant strains.xlsx', 'Foglio1');

% Get list of ORFs from both sets of ORFs
hit_orfs = data(4:end,1);
all_orfs = dataAll(2:end, 1);

% Get rid of empty rows
inds = find(cellfun(@isnumeric, all_orfs));
all_orfs(inds) = [];

% Clean up ORFs
hit_orfs = cleanGenename(hit_orfs);
all_orfs = cleanGenename(all_orfs);

% Standardize ORFs
hit_orfs = translate(hit_orfs);
all_orfs = translate(all_orfs);

% Find anything that doesn't look like an ORF
inds = find(~isorf(hit_orfs));
disp(hit_orfs(inds)); 
inds = find(~isorf(all_orfs));
disp(all_orfs(inds));

% Make sure all are unique
all_orfs = unique(all_orfs);

% Make a zero matrix for all data points
all_data = zeros(length(all_orfs), 2);

% Get data from hits and clean it up
hit_data = data(4:end, 8:9);

% Set Numeric Values to hits data
hit_data(strcmp('HS', hit_data)) = {-3};
hit_data(strcmp('MS', hit_data)) = {-2};
hit_data(strcmp('LS', hit_data)) = {-1};

% Transform cell array into Double Array
hit_data = cell2mat(hit_data);

% Average any repeated value
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% Check to see all hit_orfs are in all_orfs
% missing_orfs = setdiff(hit_orfs, all_orfs); 
% disp(missing_orfs);

% Find indices for hit_orfs in all_orfs 
[~, ind1, ind2] = intersect(all_orfs, hit_orfs);
all_data(ind1) = hit_data(ind2);

% Prepare final dataset
tigano_ottonello_2015.orfs = all_orfs;
tigano_ottonello_2015.ph = strcat(phenotypes, '; ', treatments);
tigano_ottonello_2015.data = all_data;

%% Save

save('./tigano_ottonello_2015.mat','tigano_ottonello_2015');

end
