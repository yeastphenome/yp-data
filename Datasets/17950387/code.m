%% Galvan~Smith, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
galvan_smith_2008.pmid = 17950387;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(galvan_smith_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load plate maps

[FILENAMES{end+1}, map.raw] = read_data('xlsread','./raw_data/yGDA-Master_Plate_list_Combined(New).xlsx', 'Sheet1');
map.raw(1,:) = [];
map.platerowcol = cell2mat(map.raw(:,4:6));

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/IMELDA 08Feb2006GDAraw data.xls', 'Sheet1');
data.raw(1,:) = [];
data.platerowcol = cell2mat(data.raw(:,1:3));
data.orfs = cell(size(data.raw,1),1);

[C,ia,ib] = intersect(data.platerowcol,map.platerowcol,'rows');
data.orfs(ia) = map.raw(ib,2);

data.scores = cell2mat(data.raw(:,10:11));
data.scores_norm = data.scores(:,2)./data.scores(:,1);  % Test normalized divided by control normalized

data.orfs = clean_orf(data.orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(data.orfs));
disp(data.orfs(inds)); 

data.orfs(inds) = [];
data.scores_norm(inds) = [];

% Average data for identical ORFs that appear multiple times
[data.orfs, data.scores_norm] = grpstats(data.scores_norm, data.orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [134];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
galvan_smith_2008.orfs = data.orfs;
galvan_smith_2008.ph = hit_data_names;
galvan_smith_2008.data = data.scores_norm;
galvan_smith_2008.dataset_ids = hit_data_ids;

%% Save

save('./galvan_smith_2008.mat','galvan_smith_2008');

%% Print out

fid = fopen('./galvan_smith_2008.txt','w');
write_matrix_file(fid, galvan_smith_2008.orfs, galvan_smith_2008.ph, galvan_smith_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(galvan_smith_2008)
end

end
