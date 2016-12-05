%% Okada~Shima, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
okada_shima_2014.pmid = 24410772;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(okada_shima_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,3);
hit_data = cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [570];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
okada_shima_2014.orfs = hit_strains;
okada_shima_2014.ph = hit_data_names;
okada_shima_2014.data = hit_data;
okada_shima_2014.dataset_ids = hit_data_ids;

%% Save

save('./okada_shima_2014.mat','okada_shima_2014');

%% Print out

fid = fopen('./okada_shima_2014.txt','w');
write_matrix_file(fid, okada_shima_2014.orfs, okada_shima_2014.ph, okada_shima_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(okada_shima_2014)
end

end

