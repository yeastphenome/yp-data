%% Yadav~Yadav, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yadav_yadav_2011.pmid = 20967895;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yadav_yadav_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Workbook1yea_1825_supportinginfor.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains_1 = data(2,1);
hit_strains_2 = data(2,2);
hit_strains_3 = data(2,3);

hit_strains_1 = strsplit(hit_strains_1{1}, ',');
hit_strains_2 = strsplit(hit_strains_2{1}, ',');
hit_strains_3 = strsplit(hit_strains_3{1}, ',');

hit_strains_1 = hit_strains_1';
hit_strains_1 = clean_genename(hit_strains_1);
hit_strains_1 = translate(hit_strains_1);
hit_strains_2 = hit_strains_2';
hit_strains_2 = clean_genename(hit_strains_2);
hit_strains_2 = translate(hit_strains_2);
hit_strains_3 = hit_strains_3';
hit_strains_3 = clean_genename(hit_strains_3);
hit_strains_3 = translate(hit_strains_3);

hit_strains = [hit_strains_1; hit_strains_2; hit_strains_3];
hit_strains = unique(hit_strains);

% Get the data itself
hit_data = nan(length(hit_strains), 3);
[~, ind1, ind2] = intersect(hit_strains, hit_strains_1);
hit_data(ind1, 1) = -1;
[~, ind1, ind2] = intersect(hit_strains, hit_strains_2);
hit_data(ind1, 2) = -1;
[~, ind1, ind2] = intersect(hit_strains, hit_strains_3);
hit_data(ind1, 3) = -1;

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds, :) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [132; 446; 447];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yadav_yadav_2011.orfs = hit_strains;
yadav_yadav_2011.ph = hit_data_names;
yadav_yadav_2011.data = hit_data;
yadav_yadav_2011.dataset_ids = hit_data_ids;

%% Save

save('./yadav_yadav_2011.mat','yadav_yadav_2011');

%% Print out

fid = fopen('./yadav_yadav_2011.txt','w');
write_matrix_file(fid, yadav_yadav_2011.orfs, yadav_yadav_2011.ph, yadav_yadav_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yadav_yadav_2011)
end

end

