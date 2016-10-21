%% Svensson~Samson, 2013
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
svensson_samson_2013.pmid = 24040048;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(svensson_samson_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/logratios(treatedcontrol)-3807strains.xlsx');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
strains(inds) = [];

%% Get data from hits
% Make a data matrix
hit_data = data(2:end, 3:end);
hit_data = cell2mat(hit_data);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5217; 5221; 5220; 5219; 5222; 5223; 5224; 5218; 5225; 5229; 5228; 5227; 5230; 5231; 5232; 5226];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

svensson_samson_2013.orfs = strains;
svensson_samson_2013.data = hit_data;
svensson_samson_2013.ph = hit_data_names;
svensson_samson_2013.dataset_ids = hit_data_ids;

%% Save

save('./svensson_samson_2013.mat','svensson_samson_2013');

%% Print out

fid = fopen('./svensson_samson_2013.txt','w');
write_matrix_file(fid, svensson_samson_2013.orfs, svensson_samson_2013.ph, svensson_samson_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(svensson_samson_2013)
end

end
