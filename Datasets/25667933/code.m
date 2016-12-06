%% Nislow~Hammond, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
nislow_hammond_2015.pmid = 25667933;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(nislow_hammond_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data - first file: ground

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS5_vsT1-hom1-dropsToBg.xlsx', 'ground');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_1 = data(2:end,2);

% Get the data itself
hit_data_1 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_1 = clean_orf(hit_strains_1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_1));
hit_strains_1(inds) = [];
hit_data_1(inds, :) = [];

hit_data_1 = cell2mat(hit_data_1);

% If the same strain is present more than once, average its values
[hit_strains_1, hit_data_1] = grpstats(hit_data_1, hit_strains_1, {'gname','mean'});

%% Load the data - first file: flight

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS5_vsT1-hom1-dropsToBg.xlsx', 'flight');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_2 = data(2:end,2);

% Get the data itself
hit_data_2 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_2 = clean_orf(hit_strains_2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_2));
hit_strains_2(inds) = [];
hit_data_2(inds, :) = [];

hit_data_2 = cell2mat(hit_data_2);

% If the same strain is present more than once, average its values
[hit_strains_2, hit_data_2] = grpstats(hit_data_2, hit_strains_2, {'gname','mean'});

%% Load the data - second file: ground

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS6_vsT1-hom1-NaCl-dropsToBg.xlsx', 'ground');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_3 = data(2:end,2);

% Get the data itself
hit_data_3 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_3 = clean_orf(hit_strains_3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_3));
hit_strains_3(inds) = [];
hit_data_3(inds, :) = [];

hit_data_3 = cell2mat(hit_data_3);

% If the same strain is present more than once, average its values
[hit_strains_3, hit_data_3] = grpstats(hit_data_3, hit_strains_3, {'gname','mean'});

%% Load the data - second file: flight

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS6_vsT1-hom1-NaCl-dropsToBg.xlsx', 'flight');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_4 = data(2:end,2);

% Get the data itself
hit_data_4 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_4 = clean_orf(hit_strains_4);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_4));
hit_strains_4(inds) = [];
hit_data_4(inds, :) = [];

hit_data_4 = cell2mat(hit_data_4);

% If the same strain is present more than once, average its values
[hit_strains_4, hit_data_4] = grpstats(hit_data_4, hit_strains_4, {'gname','mean'});

%% Load the data - third file: ground

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS8_vsT1-het1-dropsToBg.xlsx', 'ground');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_5 = data(2:end,2);

% Get the data itself
hit_data_5 = data(2:end,4);
   
% Eliminate all white spaces & capitalize
hit_strains_5 = clean_orf(hit_strains_5);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_5));
hit_strains_5(inds) = [];
hit_data_5(inds, :) = [];

hit_data_5 = cell2mat(hit_data_5);

% If the same strain is present more than once, average its values
[hit_strains_5, hit_data_5] = grpstats(hit_data_5, hit_strains_5, {'gname','mean'});

%% Load the data - third file: flight

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS8_vsT1-het1-dropsToBg.xlsx', 'flight');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_6 = data(2:end,2);

% Get the data itself
hit_data_6 = data(2:end,4);
   
% Eliminate all white spaces & capitalize
hit_strains_6 = clean_orf(hit_strains_6);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_6));
hit_strains_6(inds) = [];
hit_data_6(inds, :) = [];

hit_data_6 = cell2mat(hit_data_6);

% If the same strain is present more than once, average its values
[hit_strains_6, hit_data_6] = grpstats(hit_data_6, hit_strains_6, {'gname','mean'});

%% Load the data - fourth file: ground

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS9_vsT1-het1-NaCl-dropsToBg.xlsx', 'ground');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_7 = data(2:end,2);

% Get the data itself
hit_data_7 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_7 = clean_orf(hit_strains_7);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_7));
hit_strains_7(inds) = [];
hit_data_7(inds, :) = [];

hit_data_7 = cell2mat(hit_data_7);

% If the same strain is present more than once, average its values
[hit_strains_7, hit_data_7] = grpstats(hit_data_7, hit_strains_7, {'gname','mean'});

%% Load the data - fourth file: flight

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TableS9_vsT1-het1-NaCl-dropsToBg.xlsx', 'flight');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_8 = data(2:end,2);

% Get the data itself
hit_data_8 = data(2:end,4:5);
   
% Eliminate all white spaces & capitalize
hit_strains_8 = clean_orf(hit_strains_8);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_8));
hit_strains_8(inds) = [];
hit_data_8(inds, :) = [];

hit_data_8 = cell2mat(hit_data_8);

% If the same strain is present more than once, average its values
[hit_strains_8, hit_data_8] = grpstats(hit_data_8, hit_strains_8, {'gname','mean'});

%% Combine the data

hit_strains = [hit_strains_1; hit_strains_2; hit_strains_3; hit_strains_4; hit_strains_5; hit_strains_6; hit_strains_7; hit_strains_8];
hit_strains = unique(hit_strains);

hit_data = nan(length(hit_strains), 14);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_1);
hit_data(ind1, 1) = hit_data_1(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_1);
hit_data(ind1, 2) = hit_data_1(ind2,2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_2);
hit_data(ind1, 3) = hit_data_2(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_2);
hit_data(ind1, 4) = hit_data_2(ind2,2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_3);
hit_data(ind1, 5) = hit_data_3(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_3);
hit_data(ind1, 6) = hit_data_3(ind2,2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_4);
hit_data(ind1, 7) = hit_data_4(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_4);
hit_data(ind1, 8) = hit_data_4(ind2,2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_5);
hit_data(ind1, 9) = hit_data_5(ind2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_6);
hit_data(ind1, 10) = hit_data_6(ind2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_7);
hit_data(ind1, 11) = hit_data_7(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_7);
hit_data(ind1, 12) = hit_data_7(ind2,2);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_8);
hit_data(ind1, 13) = hit_data_8(ind2,1);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_8);
hit_data(ind1, 14) = hit_data_8(ind2,2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5276; 5277; 5278; 5279; 5280; 5281; 5282; 5283; 5284; 5285; 5286; 5287; 5288; 5289];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
nislow_hammond_2015.orfs = hit_strains;
nislow_hammond_2015.ph = hit_data_names;
nislow_hammond_2015.data = hit_data;
nislow_hammond_2015.dataset_ids = hit_data_ids;

%% Save

save('./nislow_hammond_2015.mat','nislow_hammond_2015');

%% Print out

fid = fopen('./nislow_hammond_2015.txt','w');
write_matrix_file(fid, nislow_hammond_2015.orfs, nislow_hammond_2015.ph, nislow_hammond_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(nislow_hammond_2015)
end

end

