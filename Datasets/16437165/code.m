%% Ravid~Hochstrasser, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ravid_hochstrasser_2006.pmid = 16437165;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ravid_hochstrasser_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(:,1);

% Get the data itself
hit_data = cell2mat(data(:,2));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [430];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ravid_hochstrasser_2006.orfs = hit_strains;
ravid_hochstrasser_2006.ph = hit_data_names;
ravid_hochstrasser_2006.data = hit_data;
ravid_hochstrasser_2006.dataset_ids = hit_data_ids;

%% Save

save('./ravid_hochstrasser_2006.mat','ravid_hochstrasser_2006');

%% Print out

fid = fopen('./ravid_hochstrasser_2006.txt','w');
write_matrix_file(fid, ravid_hochstrasser_2006.orfs, ravid_hochstrasser_2006.ph, ravid_hochstrasser_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ravid_hochstrasser_2006)
end

end

