%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
lum_shoemaker_2004.pmid = 14718172;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lum_shoemaker_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'P values');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);
hit_datasets = data(1,3:end)';

% Get the data itself
hit_data = data(2:end,3:end); % if the dataset is discrete or binary
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

inds = find(strcmp('NaN', hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[FILENAMES{end+1}, data] = read_data('textread','./extras/dataset_ids.txt', '%d %s','delimiter','\t');
[~,ind1,ind2] = intersect(data{2}, hit_datasets);
hit_data_ids(ind2) = data{1}(ind1);
hit_data_ids = hit_data_ids';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lum_shoemaker_2004.orfs = hit_strains;
lum_shoemaker_2004.ph = hit_data_names;
lum_shoemaker_2004.data = hit_data;
lum_shoemaker_2004.dataset_ids = hit_data_ids;

%% Save

save('./lum_shoemaker_2004.mat','lum_shoemaker_2004');

%% Print out

fid = fopen('./lum_shoemaker_2004.txt','w');
write_matrix_file(fid, lum_shoemaker_2004.orfs, lum_shoemaker_2004.ph, lum_shoemaker_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lum_shoemaker_2004)
end

end

