%% Askree~McEachern, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
askree_mceachern_2004.pmid = 15161972;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(askree_mceachern_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/01263Table3.xlsx', 'Sheet1');

% Get indices of the data columns
ind_data = 3:4;

hit_orfs = data.raw(:,1);
hit_data = data.raw(:,ind_data);

% Set data for inconsistent mutants to NaN
inds = find(strncmp('*Y', hit_orfs,2));
hit_data(inds,:) = {NaN};

% Eliminate anything that doesn't look like an ORF
hit_orfs = clean_orf(hit_orfs);
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

hit_orfs(inds) = [];
hit_data(inds,:) = [];

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};

hit_data = nanmean(cell2mat(hit_data),2);

% Average data for identical ORFs that appear multiple times
[hit_orfs,hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [113];

%% Load the tested genes

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/S.cerftpmata.xlsx', 'Sheet1');

tested_orfs = data.raw(:,2);
tested_orfs = clean_orf(tested_orfs);

tested_orfs(strcmp('YOLO57W', tested_orfs)) = {'YOL057W'};
tested_orfs(strcmp('YOLO62C', tested_orfs)) = {'YOL062C'};
tested_orfs(strcmp('YKLO72W', tested_orfs)) = {'YKL072W'};
tested_orfs(strcmp('YJL206-A', tested_orfs)) = {'YJL206C-A'};
tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};

inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

askree_mceachern_2004.orfs = tested_orfs;
askree_mceachern_2004.ph = hit_data_names;
askree_mceachern_2004.data = zeros(length(tested_orfs),length(hit_data_names));
[~,ind1,ind2] = intersect(hit_orfs, askree_mceachern_2004.orfs);
askree_mceachern_2004.data(ind2,:) = hit_data(ind1,:);
askree_mceachern_2004.dataset_ids = hit_data_ids;

%% Save

save('./askree_mceachern_2004.mat','askree_mceachern_2004');

%% Print out

fid = fopen('./askree_mceachern_2004.txt','w');
write_matrix_file(fid, askree_mceachern_2004.orfs, askree_mceachern_2004.ph, askree_mceachern_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(askree_mceachern_2004)
end

end
