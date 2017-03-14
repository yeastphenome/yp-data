%% Vahey~Voldman, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
vahey_voldman_2013.pmid = 23661198;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(vahey_voldman_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load Hit Strains
% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/c3lc50162k.xlsx', 'Sheet1');

% Get the list of ORFs
strains = data(:,1);
hit_data = data(:,3:4);

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

strains(inds,:) = [];
hit_data(inds,:) = [];

% Get data from hits
hit_data = cell2mat(hit_data);

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

hit_data = nanmean(hit_data,2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [95];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
vahey_voldman_2013.orfs = strains;
vahey_voldman_2013.ph = hit_data_names;
vahey_voldman_2013.data = hit_data;
vahey_voldman_2013.dataset_ids = hit_data_ids;

%% Save

save('./vahey_voldman_2013.mat','vahey_voldman_2013');

%% Print out

fid = fopen('./vahey_voldman_2013.txt','w');
write_matrix_file(fid, vahey_voldman_2013.orfs, vahey_voldman_2013.ph, vahey_voldman_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(vahey_voldman_2013)
end

end
