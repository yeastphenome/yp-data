%% Huang~O'Shea, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

huang_oshea_2005.pmid = 15695358;

phenotypes = {'Pho5 phosphatase activity'};
treatments = {'phosphate starvation'};

plates = [1:51,70:71];
rows = 'ABCDEFGH';
cols = 1:12;
timepoints = [0:120:360];

[FILENAMES{end+1}, map.raw] = readdata('xlsread','./raw_data/MATa Collection.xlsx', 'MATa Collection.xls');

map.orf = map.raw(5:end,2);
map.plate = cell2mat(map.raw(5:end, 5));
map.row = map.raw(5:end,6);
map.col = cell2mat(map.raw(5:end,7));

inds = find(cellfun(@isnan, map.row));
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

for i = 1 : length(map.row)
map.row{i} = findstr(map.row{i}, rows);
end

map.row = cell2mat(map.row);
map.idx = sub2ind([8 12],map.row, map.col);
map.data = zeros(length(map.orf),4);

inds = cellfun(@isnumeric, map.orf);

map.orf = cellfun(@strtrim, map.orf,'UniformOutput',0);
inds = ~strncmp('Y', map.orf,1);


[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/Initial Screen Data.xls', 'Sheet1');

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

huang_oshea_2005.orfs = map.orf_u;
huang_oshea_2005.data = map.data_avg_avg;
huang_oshea_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('./huang_oshea_2005.mat','huang_oshea_2005');
return;

% Save data into database
dt = huang_oshea_2005;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

