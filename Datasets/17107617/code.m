%% Freimoser~Amrhein, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
freimoser_amrhein_2006.pmid = 17107617;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(freimoser_amrhein_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

folder = './raw_data/';

[FILENAMES{end+1}, map.raw] = read_data('xlsread',[folder 'Strains_plates_annotation.xlsx']);

map.orf = map.raw(4:end, 4);
map.plate = map.raw(4:end,1);
map.row = map.raw(4:end,2);
map.col = map.raw(4:end,3);

inds = cellfun(@isnumeric, map.orf);
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

map.orf = clean_orf(map.orf);

inds = find(~is_orf(map.orf));
map.orf(inds) = [];
map.plate(inds) = [];
map.row(inds) = [];
map.col(inds) = [];

map.plate = cell2mat(map.plate);
map.col = cell2mat(map.col);

files = dir(folder);
file_names = {files.name};
file_names = file_names(strncmp('Plate', file_names,5));

[FILENAMES{end+1}, C] = read_data('textscan', './raw_data/columns_to_extract.txt', '%s %s %s %s %s %s','delimiter','\t');

all_data.plate = {};
all_data.row = {};
all_data.col = {};
all_data.data = {};
for i = 1 : length(file_names)
    [~, sheets] = xlsfinfo([folder file_names{i}]);
    ind_file = find(strcmp(file_names{i}, C{1}));

    for j = 1 : length(sheets)
        ind_sheet = find(strcmp(sheets{j}, C{2}(ind_file)));

        [FILENAMES{end+1}, data.raw] = read_data('xlsread',[folder file_names{i}], sheets{j});

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

% Find anything that doesn't look like an ORF
inds = find(~is_orf(all_data.orf));
disp(all_data.orf(inds));  

% Average multiple occurences of the same ORF
[all_data.orf_u, all_data.data_avg] = grpstats(all_data.data,all_data.orf,{'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [119];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
freimoser_amrhein_2006.orfs = all_data.orf_u;
freimoser_amrhein_2006.ph = hit_data_names;
freimoser_amrhein_2006.data = all_data.data_avg;
freimoser_amrhein_2006.dataset_ids = hit_data_ids;

%% Save

save('./freimoser_amrhein_2006.mat','freimoser_amrhein_2006');

%% Print out

fid = fopen('./freimoser_amrhein_2006.txt','w');
write_matrix_file(fid, freimoser_amrhein_2006.orfs, freimoser_amrhein_2006.ph, freimoser_amrhein_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(freimoser_amrhein_2006)
end

end
