%% Jiang~Zhang, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jiang_zhang_2014.pmid = 25331360;

hit_data_ids = [771];

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jiang_zhang_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% All Strains
[FILENAMES{end+1}, all_orfs] = read_data('xlsread','./raw_data/DELETION LIBRARY.xlsx', 'DELETION LIBRARY');

% Get the list of orfs
tested_orfs = all_orfs(3:end,2);

% Clean them up
tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YELOO1C', tested_orfs)) = {'YEL001C'};

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Hit Strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/fyr12220-sup-0003-TableS2.xlsx', 'Sheet1');

% Get the list of ORFs
hit_orfs = data(4:end, 1);
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];

%% Data
% Make an array of zeros
final_data = zeros(size(tested_orfs));

% Match the hit names with all the names
[~,ind1,ind2] = intersect(tested_orfs, hit_orfs);
final_data(ind1) = -1;

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

jiang_zhang_2014.orfs = tested_orfs;
jiang_zhang_2014.ph = hit_data_names;
jiang_zhang_2014.data = final_data;
jiang_zhang_2014.dataset_ids = hit_data_ids;

%% Save

save('./jiang_zhang_2014.mat','jiang_zhang_2014');

fid = fopen('./jiang_zhang_2014.txt','w');
write_matrix_file(fid, jiang_zhang_2014.orfs, jiang_zhang_2014.ph, jiang_zhang_2014.data);
fclose(fid);

end


