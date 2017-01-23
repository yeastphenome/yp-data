%% Loukin~Saimi, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
loukin_saimi_2007.pmid = 17314135;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(loukin_saimi_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/06-7898SupplementalTable2.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,3);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% Fix the data based on WT
hit_data = cell2mat(hit_data);
hit_data = hit_data / 420;
hit_data = 1 ./ hit_data;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1190];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
loukin_saimi_2007.orfs = hit_strains;
loukin_saimi_2007.ph = hit_data_names;
loukin_saimi_2007.data = hit_data;
loukin_saimi_2007.dataset_ids = hit_data_ids;

%% Save

save('./loukin_saimi_2007.mat','loukin_saimi_2007');

%% Print out

fid = fopen('./loukin_saimi_2007.txt','w');
write_matrix_file(fid, loukin_saimi_2007.orfs, loukin_saimi_2007.ph, loukin_saimi_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(loukin_saimi_2007)
end

end

