%% Gonzalez~Daran, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gonzalez_daran_2003.pmid = 23552365;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gonzalez_daran_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

% Load plate maps
[FILENAMES{end+1}, map.raw] = read_data('xlsread','./raw_data/Knockout collection map.xls', 'DATA');
map.raw(1,:) = [];
inds = find(cellfun(@isnumeric, map.raw(:,6)));
map.raw(inds,:) = [];
[t,ia,map.rowNum] = unique(map.raw(:,6));
map.platerowcol = [cell2mat(map.raw(:,5)) map.rowNum cell2mat(map.raw(:,7))];

% Load data
datafile = './raw_data/Butanol tolerance. Screening knockout collection.xlsx';
[status,sheets] = xlsfinfo(datafile);
all_data.plate = [];
all_data.rows = [];
all_data.cols = [];
all_data.norm_data = [];

for i = 1 : length(sheets)
    control_plate = zeros(8,12,2)+NaN;
    test_plate = zeros(8,12,2)+NaN;
    plateNum = str2num(sheets{i}(3:min(length(sheets{i}),4)));
    
    [FILENAMES{end+1}, data.raw] = read_data('xlsread',datafile, sheets{i}, 'C7:N14');
    control_plate(:,:,1) = cell2mat(data.raw);
    
    [FILENAMES{end+1}, data.raw] = read_data('xlsread',datafile, sheets{i}, 'Q7:AB14');
    control_plate(:,:,2) = cell2mat(data.raw);
    
    [FILENAMES{end+1}, data.raw] = read_data('xlsread',datafile, sheets{i}, 'C19:N26');
    test_plate(:,:,2) = cell2mat(data.raw);
    
    [FILENAMES{end+1}, data.raw] = read_data('xlsread',datafile, sheets{i}, 'Q19:AB26');
    test_plate(:,:,2) = cell2mat(data.raw);
    
    norm_data = nanmean(test_plate,3)./nanmean(control_plate,3);
    [cols,rows] = meshgrid(1:12,1:8);
    
    all_data.plate = [all_data.plate; zeros(length(norm_data(:)),1)+plateNum];
    all_data.rows = [all_data.rows; rows(:)];
    all_data.cols = [all_data.cols; cols(:)];
    all_data.norm_data = [all_data.norm_data; norm_data(:)];
    i;
end

all_data.platerowcol = [all_data.plate all_data.rows all_data.cols];

% Average replicates
all_data.platerowcol2 = cell(length(all_data.plate),1);
for i = 1 : length(all_data.plate)
    all_data.platerowcol2{i} = [num2str(all_data.plate(i)),'-', num2str(all_data.rows(i)),'-', num2str(all_data.cols(i))];
end
[t,t2] = grpstats(all_data.norm_data, all_data.platerowcol2, {'gname','mean'});

all_data2.norm_data = t2;

tmp = regexp(t,'-','split');
tmp2 = vertcat(tmp{:});
all_data2.platerowcol = cellfun(@str2num, tmp2);

% Not all positions that are in the DATA have corresponding entries in the
% MAP. But that is somewhat expected, because not all positions on the
% plates are filled with strains.

[C,ia,ib] = intersect(all_data2.platerowcol,map.platerowcol,'rows');
all_data2.orfs = cell(size(all_data2.platerowcol,1),1);
all_data2.orfs(ia) = map.raw(ib,2);

% Clean the orfs
all_data2.orfs = clean_orf(all_data2.orfs);

all_data2.orfs(strcmp('YLR287-A', all_data2.orfs)) = {'YLR287C-A'};

% Eliminate indices without ORFs
inds = find(~is_orf(all_data2.orfs));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(all_data2.norm_data, all_data2.orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [127];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
gonzalez_daran_2003.orfs = t;
gonzalez_daran_2003.ph = hit_data_names;
gonzalez_daran_2003.data = t2;
gonzalez_daran_2003.dataset_ids = hit_data_ids;

%% Save

save('./gonzalez_daran_2003.mat','gonzalez_daran_2003');

%% Print out

fid = fopen('./gonzalez_daran_2003.txt','w');
write_matrix_file(fid, gonzalez_daran_2003.orfs, gonzalez_daran_2003.ph, gonzalez_daran_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gonzalez_daran_2003)
end

end
