%% Huang~O'Shea, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
huang_oshea_2005.pmid = 15695358;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(huang_oshea_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

plates = [1:51,70:71];
rows = 'ABCDEFGH';
cols = 1:12;
timepoints = [0:120:360];

[FILENAMES{end+1}, map.raw] = read_data('xlsread','./raw_data/MATa Collection.xlsx', 'MATa Collection.xls');

map.orf = map.raw(5:end,2);
map.plate = cell2mat(map.raw(5:end, 5));
map.row = map.raw(5:end,6);
map.col = cell2mat(map.raw(5:end,7));

map.orf = clean_orf(map.orf);
map.orf(strcmp('YLR287-A', map.orf)) = {'YLR287C-A'};

inds = find(~is_orf(map.orf) | cellfun(@isnan, map.row));
disp(map.orf(inds));

map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

map.row_n = nan(size(map.row));
for i = 1 : length(map.row)
    map.row_n(i) = findstr(map.row{i}, rows);
end

map.row = map.row_n;
map.idx = sub2ind([8 12],map.row, map.col);
map.data = zeros(length(map.orf),4);

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Initial Screen Data.xls', 'Sheet1');

[all_cols, all_rows] = meshgrid(1:12,1:8);
all_idx = sub2ind([8 12], all_rows, all_cols)';
[~,ix] = sort(all_idx(:));

for i = 1 : length(plates)
    inds = find(map.plate == plates(i));
    for t = 1 : length(timepoints)
        this_data = cell2mat(data.raw(6+(i-1)*7+(t-1),2:end));
        map.data(inds,t) = this_data(ix(map.idx(inds)));
    end
end

map.data_avg = nanmean(map.data,2);

% Average multiple occurences of the same ORF
[map.orf_u, map.data_avg_avg] = grpstats(map.data_avg,map.orf,{'gname','mean'});

hit_data_ids = [139];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

huang_oshea_2005.orfs = map.orf_u;
huang_oshea_2005.data = map.data_avg_avg;
huang_oshea_2005.ph = hit_data_names;
huang_oshea_2005.dataset_ids = hit_data_ids;

%% Save

save('./huang_oshea_2005.mat','huang_oshea_2005');

%% Print out

fid = fopen('./huang_oshea_2005.txt','w');
write_matrix_file(fid, huang_oshea_2005.orfs, huang_oshea_2005.ph, huang_oshea_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(huang_oshea_2005)
end

end
