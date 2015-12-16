%% Gonzalez-Ramos~Daran, 2013
function FILENAMES = code()
FILENAMES = {};

gonzalez_ramos_daran_2003.pmid = 23552365;

phenotypes = {'growth'};
treatments = {'butanol'};

% Load plate maps
[FILENAMES{end+1}, map.raw] = dataread('xlsread','./raw_data/Knockout collection map.xls', 'DATA');
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

[FILENAMES{end+1}, data.raw] = dataread('xlsread',datafile, sheets{i}, 'C7:N14');
control_plate(:,:,1) = cell2mat(data.raw);

[FILENAMES{end+1}, data.raw] = dataread('xlsread',datafile, sheets{i}, 'Q7:AB14');
control_plate(:,:,2) = cell2mat(data.raw);

[FILENAMES{end+1}, data.raw] = dataread('xlsread',datafile, sheets{i}, 'C19:N26');
test_plate(:,:,2) = cell2mat(data.raw);

[FILENAMES{end+1}, data.raw] = dataread('xlsread',datafile, sheets{i}, 'Q19:AB26');
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

% Eliminate empty positions
inds = find(cellfun(@isempty, all_data2.orfs));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

inds = find(cellfun(@isnumeric, all_data2.orfs));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

all_data2.orfs = cellfun(@strtrim, all_data2.orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', all_data2.orfs,1));
all_data2.orfs(inds) = [];
all_data2.platerowcol(inds,:) = [];
all_data2.norm_data(inds) = [];

all_data2.orfs = upper(all_data2.orfs);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(all_data2.norm_data, all_data2.orfs, {'gname','mean'});
gonzalez_ramos_daran_2003.orfs = t;
gonzalez_ramos_daran_2003.data = t2;
gonzalez_ramos_daran_2003.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'gonzalez_ramos_daran_2003.mat'],'gonzalez_ramos_daran_2003');
return;

% Save data into database
dt = gonzalez_ramos_daran_2003;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

