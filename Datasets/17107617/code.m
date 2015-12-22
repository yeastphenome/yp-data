%% Freimoser~Amrhein, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

freimoser_amrhein_2006.pmid = 17107617;

phenotypes = {'inorganic polyphosphate abundance'};
treatments = {''};

folder = './raw_data/';

[FILENAMES{end+1}, map.raw] = dataread('xlsread',[folder 'Strains_plates_annotation.xlsx']);

map.orf = map.raw(4:end, 4);
map.plate = map.raw(4:end,1);
map.row = map.raw(4:end,2);
map.col = map.raw(4:end,3);

inds = cellfun(@isnumeric, map.orf);
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

inds = ~strncmp('Y', map.orf,1);
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

map.plate = cell2mat(map.plate);
map.col = cell2mat(map.col);

files = dir([folder 'Archive/']);
file_names = {files.name};
file_names = file_names(strncmp('Plate', file_names,5));

fid = fopen('./raw_data/columns_to_extract.txt');
C = textscan(fid,'%s %s %s %s %s %s','delimiter','\t');
fclose(fid);

all_data.plate = {};
all_data.row = {};
all_data.col = {};
all_data.data = {};
for i = 1 : length(file_names)
[~, sheets] = xlsfinfo([folder 'Archive/' file_names{i}]);
ind_file = find(strcmp(file_names{i}, C{1}));

for j = 1 : length(sheets)
ind_sheet = find(strcmp(sheets{j}, C{2}(ind_file)));

[FILENAMES{end+1}, data.raw] = dataread('xlsread',[folder 'Archive/' file_names{i}], sheets{j});



ind_plate = strfind('ABCDEFGHJIKLMNOP', C{3}{ind_file(ind_sheet)});
ind_row = strfind('ABCDEFGHJIKLMNOP', C{4}{ind_file(ind_sheet)});
ind_col = strfind('ABCDEFGHJIKLMNOP', C{5}{ind_file(ind_sheet)});
ind_data = strfind('ABCDEFGHJIKLMNOP', C{6}{ind_file(ind_sheet)});

if isempty(ind_plate) || isempty(ind_row) || isempty(ind_col)
fprintf('%s\t%s\n', file_names{i}, sheets{j});
else
all_data.plate = [all_data.plate; data.raw(:,ind_plate)];
all_data.row = [all_data.row; data.raw(:,ind_row)];
all_data.col = [all_data.col; data.raw(:,ind_col)];
all_data.data = [all_data.data; data.raw(:,end-2)];
end
end
end

all_data.plate = cell2mat(all_data.plate);
all_data.col = cell2mat(all_data.col);

inds = isnan(all_data.plate) | isnan(all_data.col) | cellfun(@isnumeric, all_data.row) | ~cellfun(@isnumeric, all_data.data);

all_data.plate(inds) = [];
all_data.row(inds) = [];
all_data.col(inds) = [];
all_data.data(inds) = [];

all_data.data = cell2mat(all_data.data);


% Map ORFs
all_data.orf = cell(size(all_data.data));
for i = 1 : length(all_data.plate)
inds = find(map.plate == all_data.plate(i) & map.col == all_data.col(i) & strcmp(all_data.row{i}, map.row));
if ~isempty(inds)
all_data.orf(i,1) = map.orf(inds);
end
end

inds = find(cellfun(@isempty, all_data.orf));
all_data.orf(inds) = [];
all_data.plate(inds) = [];
all_data.row(inds) = [];
all_data.col(inds) = [];
all_data.data(inds) = [];


% Average multiple occurences of the same ORF
[all_data.orf_u, all_data.data_avg] = grpstats(all_data.data,all_data.orf,{'gname','mean'});

freimoser_amrhein_2006.orfs = all_data.orf_u;
freimoser_amrhein_2006.data = all_data.data_avg;
freimoser_amrhein_2006.ph = [strcat(phenotypes, '; ', treatments)];

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'freimoser_amrhein_2006.mat'],'freimoser_amrhein_2006');
return;

% Save data into database
dt = freimoser_amrhein_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

