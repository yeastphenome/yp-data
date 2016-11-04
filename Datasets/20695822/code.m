%% Landstetter~Kuchler, 2010
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
landstetter_kuchler_2010.pmid = 20695822;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(landstetter_kuchler_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Yeast_deletions_set.xlsx', 'Stammliste');

% Get the data
data = data(2:4828, 2:end);

% Get the orfs and clean them up
strains = data(:,1);
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
strains(inds) = [];

%% Get data from hits
% Make a data matrix
[FILENAMES{end+1}, hits] = read_data('xlsread','./raw_data/Supp_Table2.xlsx');
hit_strains = hits(6:end, 2);

% Remove anything that's not an orf
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [152];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

landstetter_kuchler_2010.orfs = strains;
landstetter_kuchler_2010.ph = hit_data_names;
landstetter_kuchler_2010.data = zeros(length(landstetter_kuchler_2010.orfs),length(landstetter_kuchler_2010.ph));
landstetter_kuchler_2010.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, landstetter_kuchler_2010.orfs);
landstetter_kuchler_2010.data(ind2,:) = -1;

%% Save

save('./landstetter_kuchler_2010.mat','landstetter_kuchler_2010');

%% Print out

fid = fopen('./landstetter_kuchler_2010.txt','w');
write_matrix_file(fid, landstetter_kuchler_2010.orfs, landstetter_kuchler_2010.ph, landstetter_kuchler_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(landstetter_kuchler_2010)
end

end
