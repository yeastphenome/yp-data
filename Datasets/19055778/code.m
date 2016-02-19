%% Alamgir~Golshani, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alamgir_golshani_2008.pmid = 19055778;

phenotypes = {'growth [colony size]'};
treatments = {'paromomycin'};

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

inds = cellfun(@isnumeric, data.orf);
data.orf(inds) = [];
data.data(inds) = [];

data.orf = cellfun(@strtrim, data.orf,'UniformOutput',0);

inds = ~strncmp('Y', data.orf,1);
data.orf(inds) = [];
data.data(inds) = [];
data.orf = upper(data.orf);

[data_avg.orf, data_avg.data] = grpstats(data.data,data.orf,{'gname','mean'});

alamgir_golshani_2008.orfs = data_avg.orf;
alamgir_golshani_2008.data = data_avg.data;

alamgir_golshani_2008.ph = strcat(phenotypes, '; ', treatments);

save('./alamgir_golshani_2008.mat','alamgir_golshani_2008');
return;

% Save data into database
dt = alamgir_golshani_2008;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


fid = fopen('./alamgir_golshani_2008.txt','w');
write_matrix_file(fid, alamgir_golshani_2008.orfs, alamgir_golshani_2008.ph, alamgir_golshani_2008.data);
fclose(fid);

end

