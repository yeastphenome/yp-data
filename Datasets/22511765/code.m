%% Kim~Cunningham, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
kim_cunningham_2012.pmid = 22511765;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(kim_cunningham_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/jbc.M112.363390-1.xlsx', 'sup table 1');

% Get ORFs & data
hit_orfs = data.raw(3:end, 5);
hit_data = data.raw(3:end, 14:15);

% Eliminate anything that doesn't look like an ORF
hit_orfs = clean_orf(hit_orfs);

hit_orfs(strcmp('YLR287-A', hit_orfs)) = {'YLR287C-A'};
inds = find(~is_orf(hit_orfs));
hit_orfs(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [19; 243];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

kim_cunningham_2012.orfs = hit_orfs;
kim_cunningham_2012.ph = hit_data_names;
kim_cunningham_2012.data = hit_data;
kim_cunningham_2012.dataset_ids = hit_data_ids;

%% Save

save('./kim_cunningham_2012.mat','kim_cunningham_2012');

%% Print out

fid = fopen('./kim_cunningham_2012.txt','w');
write_matrix_file(fid, kim_cunningham_2012.orfs, kim_cunningham_2012.ph, kim_cunningham_2012.data);
fclose(fid);

end
