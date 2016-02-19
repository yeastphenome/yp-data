%% Gresham~Botstein, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

gresham_botstein_2011.pmid = 20944018;

phenotypes = {'Death rate (%/hr)'};
treatments = {'phosphate starvation';'leucine starvation'};

% Dataset #1
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TABLES4.xlsx', 'phoAbs.txt');

data.genenames = data.raw(2:end, 1);
data.data = cell2mat(data.raw(2:end, 2));
data.genenames_noannot = cell(size(data.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data.genenames)
C = regexp(data.genenames{i},'_','split');
data.genenames_noannot{i} = C{1};
end

[data.orfs, translated] = translate(data.genenames_noannot);
data.orfs(~translated) = [];
data.data(~translated,:) = [];


% Dataset #2
[FILENAMES{end+1}, data2.raw] = read_data('xlsread','./raw_data/TABLES5.xlsx', 'leuAbs.txt');

data2.genenames = data2.raw(2:end, 1);
data2.data = cell2mat(data2.raw(2:end, 2));
data2.genenames_noannot = cell(size(data2.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data2.genenames)
C = regexp(data2.genenames{i},'_','split');
data2.genenames_noannot{i} = C{1};
end

[data2.orfs, translated] = translate(data2.genenames_noannot);
data2.orfs(~translated) = [];
data2.data(~translated,:) = [];


% Average data for identical ORFs that appear multiple times
[data.orfs_u,data.data_avg] = grpstats(data.data, data.orfs, {'gname','mean'});
[data2.orfs_u,data2.data_avg] = grpstats(data2.data, data2.orfs, {'gname','mean'});

gresham_botstein_2011.orfs = unique([data.orfs_u; data2.orfs_u]);
gresham_botstein_2011.data = nan(length(gresham_botstein_2011.orfs),2);

[~,ind1,ind2] = intersect(gresham_botstein_2011.orfs, data.orfs_u);
gresham_botstein_2011.data(ind1,1) = data.data_avg(ind2);
[~,ind1,ind2] = intersect(gresham_botstein_2011.orfs, data2.orfs_u);
gresham_botstein_2011.data(ind1,2) = data2.data_avg(ind2);

gresham_botstein_2011.ph = strcat(phenotypes, {'; '}, treatments);

save('./gresham_botstein_2011.mat','gresham_botstein_2011');
return;

% Save data into database
dt = gresham_botstein_2011;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


fid = fopen('./gresham_botstein_2011.txt','w');
write_matrix_file(fid, gresham_botstein_2011.orfs, gresham_botstein_2011.ph, gresham_botstein_2011.data);
fclose(fid);

end

