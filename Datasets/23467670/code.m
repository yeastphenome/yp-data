%% Samanfar~Golshani, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
samanfar_golshani_2013.pmid = 23467670;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(samanfar_golshani_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/c3mb25516f-2.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(:,2);

% Get the data itself
hit_data = data(:,3);
hit_data = cell2mat(hit_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [128];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
samanfar_golshani_2013.orfs = hit_strains;
samanfar_golshani_2013.ph = hit_data_names;
samanfar_golshani_2013.data = hit_data;
samanfar_golshani_2013.dataset_ids = hit_data_ids;

%% Save

save('./samanfar_golshani_2013.mat','samanfar_golshani_2013');

%% Print out

fid = fopen('./samanfar_golshani_2013.txt','w');
write_matrix_file(fid, samanfar_golshani_2013.orfs, samanfar_golshani_2013.ph, samanfar_golshani_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(samanfar_golshani_2013)
end

end
