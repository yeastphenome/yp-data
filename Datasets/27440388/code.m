%% Zhao~Liu, 2016
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhao_liu_2016.pmid = 27440388;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zhao_liu_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Beidong Percentage of cells with synphilin-1 aggregates-SGA-0803.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data(5:end,1);

% Get the data itself
hit_data = data(5:end,4);
indx = find(strcmp('NS', hit_data));
hit_data(indx) = {0};
indx = find(strcmp('ND', hit_data));
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

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
hit_data_ids = [4825];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zhao_liu_2016.orfs = hit_strains;
zhao_liu_2016.ph = hit_data_names;
zhao_liu_2016.data = hit_data;
zhao_liu_2016.dataset_ids = hit_data_ids;

%% Save
save('./zhao_liu_2016.mat','zhao_liu_2016');

%% Print out
fid = fopen('./zhao_liu_2016.txt','w');
write_matrix_file(fid, zhao_liu_2016.orfs, zhao_liu_2016.ph, zhao_liu_2016.data);
fclose(fid);

%% Save to DB (admin)
addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zhao_liu_2016)
end

end
