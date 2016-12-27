%% Shima~Takagi, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
shima_takagi_2008.pmid = 18224659;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(shima_takagi_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Air-drying stress.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:end,1);

% Get the data itself
hit_data = data(4:end,4);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds, :) = [];

% Remove data values that are ND(Not determined)
hit_data(strcmp(hit_data, 'ND')) = {NaN};
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [504];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
shima_takagi_2008.orfs = hit_strains;
shima_takagi_2008.ph = hit_data_names;
shima_takagi_2008.data = hit_data;
shima_takagi_2008.dataset_ids = hit_data_ids;

%% Save

save('./shima_takagi_2008.mat','shima_takagi_2008');

%% Print out

fid = fopen('./shima_takagi_2008.txt','w');
write_matrix_file(fid, shima_takagi_2008.orfs, shima_takagi_2008.ph, shima_takagi_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(shima_takagi_2008)
end

end

