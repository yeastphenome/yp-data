%% Thorsen~Tamas, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
thorsen_tamas_2009.pmid = 19284616;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(thorsen_tamas_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
% As (III) - 0.5 mM & 24 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '0.5 As 24');

% Get the list of ORFs and the correponding data 
hit_strains_1 = data(2:end,11);

% Get the data itself
hit_data_1 = data(2:end,10);
hit_data_1 = cell2mat(hit_data_1);
   
% Eliminate all white spaces & capitalize
hit_strains_1 = clean_orf(hit_strains_1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_1));
hit_strains_1(inds) = [];
hit_data_1(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_1, hit_data_1] = grpstats(hit_data_1, hit_strains_1, {'gname','mean'});

% As (III) - 1 mM & 24 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '1 As 24h');

% Get the list of ORFs and the correponding data 
hit_strains_2 = data(2:end,11);

% Get the data itself
hit_data_2 = data(2:end,10);
hit_data_2 = cell2mat(hit_data_2);
   
% Eliminate all white spaces & capitalize
hit_strains_2 = clean_orf(hit_strains_2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_2));
hit_strains_2(inds) = [];
hit_data_2(inds) = [];

% If the same strain is present more than once, average its value
[hit_strains_2, hit_data_2] = grpstats(hit_data_2, hit_strains_2, {'gname','mean'});

% As (III) - 1.5 mM & 24 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '1.5 As 24h');

% Get the list of ORFs and the correponding data 
hit_strains_3 = data(2:end,11);

% Get the data itself
hit_data_3 = data(2:end,10);
hit_data_3 = cell2mat(hit_data_3);
   
% Eliminate all white spaces & capitalize
hit_strains_3 = clean_orf(hit_strains_3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_3));
hit_strains_3(inds) = [];
hit_data_3(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_3, hit_data_3] = grpstats(hit_data_3, hit_strains_3, {'gname','mean'});

% As (III) - 0.5 mM & 48 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '0.5As 48h');

% Get the list of ORFs and the correponding data 
hit_strains_4 = data(2:end,11);

% Get the data itself
hit_data_4 = data(2:end,10); 
hit_data_4 = cell2mat(hit_data_4);
   
% Eliminate all white spaces & capitalize
hit_strains_4 = clean_orf(hit_strains_4);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_4));
hit_strains_4(inds) = [];
hit_data_4(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_4, hit_data_4] = grpstats(hit_data_4, hit_strains_4, {'gname','mean'});

% As (III) - 1 mM & 48 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '1 As 48h');

% Get the list of ORFs and the correponding data 
hit_strains_5 = data(2:end,11);

% Get the data itself
hit_data_5 = data(2:end,10);
hit_data_5 = cell2mat(hit_data_5);
   
% Eliminate all white spaces & capitalize
hit_strains_5 = clean_orf(hit_strains_5);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_5));
hit_strains_5(inds) = [];
hit_data_5(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_5, hit_data_5] = grpstats(hit_data_5, hit_strains_5, {'gname','mean'});

% As (III) - 1.5 mM & 48 hr

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/As(III) list avg.xlsx', '1.5 As 48h');

% Get the list of ORFs and the correponding data 
hit_strains_6 = data(2:end,11);

% Get the data itself
hit_data_6 = data(2:end,10); 
hit_data_6 = cell2mat(hit_data_6);
   
% Eliminate all white spaces & capitalize
hit_strains_6 = clean_orf(hit_strains_6);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_6));
hit_strains_6(inds) = [];
hit_data_6(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_6, hit_data_6] = grpstats(hit_data_6, hit_strains_6, {'gname','mean'});

% Cd - 75 uM & 24 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 75 24h');

% Get the list of ORFs and the correponding data 
hit_strains_7 = data(2:end,10);

% Get the data itself
hit_data_7 = data(2:end,4);
hit_data_7(~cellfun(@(x) any(isnumeric(x(:))), hit_data_7)) = {NaN};
hit_data_7 = cell2mat(hit_data_7);

% Eliminate all white spaces & capitalize
hit_strains_7 = clean_orf(hit_strains_7);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_7));
hit_strains_7(inds) = [];
hit_data_7(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_7, hit_data_7] = grpstats(hit_data_7, hit_strains_7, {'gname','mean'});

% Cd - 100 uM & 24 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 100 24h');

% Get the list of ORFs and the correponding data 
hit_strains_8 = data(2:end,10);

% Get the data itself
hit_data_8 = data(2:end,4);
hit_data_8(~cellfun(@(x) any(isnumeric(x(:))), hit_data_8)) = {NaN};
hit_data_8 = cell2mat(hit_data_8);
   
% Eliminate all white spaces & capitalize
hit_strains_8 = clean_orf(hit_strains_8);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_8));
hit_strains_8(inds) = [];
hit_data_8(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_8, hit_data_8] = grpstats(hit_data_8, hit_strains_8, {'gname','mean'});

% Cd - 150 uM & 24 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 150 24h');

% Get the list of ORFs and the correponding data 
hit_strains_9 = data(2:end,10);

% Get the data itself
hit_data_9 = data(2:end,4);
hit_data_9(~cellfun(@(x) any(isnumeric(x(:))), hit_data_9)) = {NaN};
hit_data_9 = cell2mat(hit_data_9);
   
% Eliminate all white spaces & capitalize
hit_strains_9 = clean_orf(hit_strains_9);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_9));
hit_strains_9(inds) = [];
hit_data_9(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_9, hit_data_9] = grpstats(hit_data_9, hit_strains_9, {'gname','mean'});

% Cd - 75 uM & 48 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 75 48h');

% Get the list of ORFs and the correponding data 
hit_strains_10 = data(2:end,10);

% Get the data itself
hit_data_10 = data(2:end,4); 
hit_data_10(~cellfun(@(x) any(isnumeric(x(:))), hit_data_10)) = {NaN};
hit_data_10 = cell2mat(hit_data_10);

% Eliminate all white spaces & capitalize
hit_strains_10 = clean_orf(hit_strains_10);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_10));
hit_strains_10(inds) = [];
hit_data_10(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_10, hit_data_10] = grpstats(hit_data_10, hit_strains_10, {'gname','mean'});

% Cd - 100 uM & 48 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 100 48h');

% Get the list of ORFs and the correponding data 
hit_strains_11 = data(2:end,10);

% Get the data itself
hit_data_11 = data(2:end,4); 
hit_data_11(~cellfun(@(x) any(isnumeric(x(:))), hit_data_11)) = {NaN};
hit_data_11 = cell2mat(hit_data_11);

% Eliminate all white spaces & capitalize
hit_strains_11 = clean_orf(hit_strains_11);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_11));
hit_strains_11(inds) = [];
hit_data_11(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_11, hit_data_11] = grpstats(hit_data_11, hit_strains_11, {'gname','mean'});

% Cd - 150 uM & 48 hr
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Cd list avg.xlsx', 'Cd 150 48h');

% Get the list of ORFs and the correponding data 
hit_strains_12 = data(2:end,10);

% Get the data itself
hit_data_12 = data(2:end,4); 
hit_data_12(~cellfun(@(x) any(isnumeric(x(:))), hit_data_12)) = {NaN};
hit_data_12 = cell2mat(hit_data_12);

% Eliminate all white spaces & capitalize
hit_strains_12 = clean_orf(hit_strains_12);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_12));
hit_strains_12(inds) = [];
hit_data_12(inds) = [];

% If the same strain is present more than once, average its values
[hit_strains_12, hit_data_12] = grpstats(hit_data_12, hit_strains_12, {'gname','mean'});

%% Combine all datasets
hit_strains = unique([hit_strains_1; hit_strains_2; hit_strains_3; hit_strains_4; hit_strains_5; hit_strains_6; hit_strains_7; hit_strains_8; hit_strains_9; hit_strains_10; hit_strains_11; hit_strains_12]);
hit_data = nan(length(hit_strains), 12);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_1);
hit_data(ind1, 1) = hit_data_1(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_2);
hit_data(ind1, 2) = hit_data_2(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_3);
hit_data(ind1, 3) = hit_data_3(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_4);
hit_data(ind1, 4) = hit_data_4(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_5);
hit_data(ind1, 5) = hit_data_5(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_6);
hit_data(ind1, 6) = hit_data_6(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_7);
hit_data(ind1, 7) = hit_data_7(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_8);
hit_data(ind1, 8) = hit_data_8(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_9);
hit_data(ind1, 9) = hit_data_9(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_10);
hit_data(ind1, 10) = hit_data_10(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_11);
hit_data(ind1, 11) = hit_data_11(ind2);

[~, ind1, ind2] = intersect(hit_strains, hit_strains_12);
hit_data(ind1, 12) = hit_data_12(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1320; 5362; 5363; 5364; 5365; 5366; 1319; 5367; 5368; 5369; 5370; 5371];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
thorsen_tamas_2009.orfs = hit_strains;
thorsen_tamas_2009.ph = hit_data_names;
thorsen_tamas_2009.data = hit_data;
thorsen_tamas_2009.dataset_ids = hit_data_ids;

%% Save

save('./thorsen_tamas_2009.mat','thorsen_tamas_2009');

%% Print out

fid = fopen('./thorsen_tamas_2009.txt','w');
write_matrix_file(fid, thorsen_tamas_2009.orfs, thorsen_tamas_2009.ph, thorsen_tamas_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(thorsen_tamas_2009)
end

end