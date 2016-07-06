%% Alamgir~Golshani, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alamgir_golshani_2008.pmid = 19055778;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(alamgir_golshani_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

% Load tested
[FILENAMES{end+1}, map.raw] = read_data('xlsread','./raw_data/Master plate list -- part 1.xlsx', 'Sheet2-Corrected');
map.orf = map.raw(2:end,2);
map.plate = cell2mat(map.raw(2:end,4));
map.row = cell2mat(map.raw(2:end,5));
map.col = cell2mat(map.raw(2:end,6));


[FILENAMES{end+1}, map2.raw] = read_data('xlsread','./raw_data/Master plate list -- part2.xlsx', 'Sheet2-Corrected');
[rinds,cinds] = find(~cellfun(@isnumeric, map2.raw(:,4:6)));
map2.raw(unique(rinds),:) = [];
map2.orf = map2.raw(:,2);
map2.plate = cell2mat(map2.raw(:,4));
map2.row = cell2mat(map2.raw(:,5));
map2.col = cell2mat(map2.raw(:,6));

map.orf = [map.orf; map2.orf];
map.plate = [map.plate; map2.plate];
map.row = [map.row; map2.row];
map.col = [map.col; map2.col];

map.inds = sub2ind([length(unique(map.plate)), length(unique(map.row)), length(unique(map.col))], map.plate, map.row, map.col);

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Plate_Analyzer_Alamgir_13mgmLParomomycin_08Dec05(1).xlsx');
data.plate = cell2mat(data.raw(2:end,1));
data.row = cell2mat(data.raw(2:end,2));
data.col = cell2mat(data.raw(2:end,3));
data.data = cell2mat(data.raw(2:end,4));

data.inds = sub2ind([length(unique(data.plate)), length(unique(data.row)), length(unique(data.col))], data.plate, data.row, data.col);
data.orf = cell(length(data.inds),1);

[~,ind1,ind2] = intersect(data.inds, map.inds);
data.orf(ind1) = map.orf(ind2);

inds = find(cellfun(@isempty, data.orf));       % 26 positions present in data don't have an ORF
data.orf(inds) = [];
data.data(inds) = [];

data.orf = clean_orf(data.orf);

data.orf(strcmp('YPL072WA', data.orf)) = {'YPL072W-A'};
inds = find(~is_orf(data.orf));
disp(data.orf(inds));

data.orf(inds) = [];
data.data(inds) = [];

[data.orf, data.data] = grpstats(data.data,data.orf,{'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [142];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

alamgir_golshani_2008.orfs = data.orf;
alamgir_golshani_2008.data = data.data;
alamgir_golshani_2008.ph = hit_data_names;
alamgir_golshani_2008.dataset_ids = hit_data_ids;

%% Save

save('./alamgir_golshani_2008.mat','alamgir_golshani_2008');

%% Print out

fid = fopen('./alamgir_golshani_2008.txt','w');
write_matrix_file(fid, alamgir_golshani_2008.orfs, alamgir_golshani_2008.ph, alamgir_golshani_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(alamgir_golshani_2008)
end

end
