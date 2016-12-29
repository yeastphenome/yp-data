%% McLaughlin~Tumer, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mclaughlin_tumer_2009.pmid = 20007368;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mclaughlin_tumer_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/0909777106_0909777106S.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,3:2:end-1);
hit_data = cell2mat(hit_data);

% Normalize by the WT
hit_data = hit_data ./ repmat(hit_data(1,:),size(hit_data,1),1) - 1;
hit_strains(1) = [];
hit_data(1,:) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds, :) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [603; 5359; 5360; 5361];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
mclaughlin_tumer_2009.orfs = hit_strains;
mclaughlin_tumer_2009.ph = hit_data_names;
mclaughlin_tumer_2009.data = hit_data;
mclaughlin_tumer_2009.dataset_ids = hit_data_ids;

%% Save

save('./mclaughlin_tumer_2009.mat','mclaughlin_tumer_2009');

%% Print out

fid = fopen('./mclaughlin_tumer_2009.txt','w');
write_matrix_file(fid, mclaughlin_tumer_2009.orfs, mclaughlin_tumer_2009.ph, mclaughlin_tumer_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mclaughlin_tumer_2009)
end

end

