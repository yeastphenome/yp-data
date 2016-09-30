%% Bin-Umer~Tumer, 2014

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bin_umer_tumer_2014.pmid = 25071194;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bin_umer_tumer_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hit Data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pnas.1403145111.st01.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,2);

% Get the data itself
hit_data = -ones(size(hit_strains));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1326];


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
bin_umer_tumer_2014.orfs = hit_strains;
bin_umer_tumer_2014.ph = hit_data_names;
bin_umer_tumer_2014.data = hit_data;
bin_umer_tumer_2014.dataset_ids = hit_data_ids;

%% Save

save('./bin_umer_tumer_2014.mat','bin_umer_tumer_2014');

%% Print out

fid = fopen('./bin_umer_tumer_2014.txt','w');
write_matrix_file(fid, bin_umer_tumer_2014.orfs, bin_umer_tumer_2014.ph, bin_umer_tumer_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bin_umer_tumer_2014)
end

end

