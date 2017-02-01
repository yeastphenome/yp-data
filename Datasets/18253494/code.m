%% Cook~Ho, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
cook_ho_2008.pmid = 18253494;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(cook_ho_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pone.0001546.s015.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,13);
hit_data = -cell2mat(hit_data); % seems that, originally, negative values = large cells, positive values = small cells. By convention, we revert that.
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5374];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
cook_ho_2008.orfs = hit_strains;
cook_ho_2008.ph = hit_data_names;
cook_ho_2008.data = hit_data;
cook_ho_2008.dataset_ids = hit_data_ids;

%% Save

save('./cook_ho_2008.mat','cook_ho_2008');

%% Print out

fid = fopen('./cook_ho_2008.txt','w');
write_matrix_file(fid, cook_ho_2008.orfs, cook_ho_2008.ph, cook_ho_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(cook_ho_2008)
end

end

