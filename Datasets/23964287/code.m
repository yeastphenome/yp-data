%% Gaytan~Vulpe, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gaytan_vulpe_2013.pmid = 23964287;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gaytan_vulpe_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/DataSheet2.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(11:end,1);

% Get the data itself
hit_data = cell2mat(data(11:end,3));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [196];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
gaytan_vulpe_2013.orfs = hit_strains;
gaytan_vulpe_2013.ph = hit_data_names;
gaytan_vulpe_2013.data = hit_data;
gaytan_vulpe_2013.dataset_ids = hit_data_ids;

%% Save

save('./gaytan_vulpe_2013.mat','gaytan_vulpe_2013');

%% Print out

fid = fopen('./gaytan_vulpe_2013.txt','w');
write_matrix_file(fid, gaytan_vulpe_2013.orfs, gaytan_vulpe_2013.ph, gaytan_vulpe_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gaytan_vulpe_2013)
end

end