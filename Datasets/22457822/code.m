%% Chesi~Gitler, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

chesi_gitler_2012.pmid = 22457822;

phenotypes = {'Growth, colony size'};
treatments = {'Manganese'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/yeast deletions Mn.xlsx', 'single deletion');

crr = zeros(size(data.raw,2)-1,3);  % Concentration, within-round replicate, round
for i = 2 : length(data.raw(1,:))
t = textscan(data.raw{1,i},'%d %s size-%d %d');
crr(i-1,:) = cell2mat(t([1 3 4]));
end

% Get indices of the data columns
ind_data = 2:19;


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = find(~strncmp('Y', data.raw(:,1),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

ind_data_conc1 = [1 2 7 8 13 14];
ind_data_conc2 = [3 4 9 10 15 16];
ind_data_conc3 = [5 6 11 12 17 18];

data2.orfs = upper(data.raw(:,1));
data2.data(:,1) = nanmean(cell2mat(t(:,ind_data_conc1)),2);
data2.data(:,2) = nanmean(cell2mat(t(:,ind_data_conc2)),2);
data2.data(:,3) = nanmean(cell2mat(t(:,ind_data_conc3)),2);

% Normalize to no drug
data2.data(:,2) = data2.data(:,2) ./ data2.data(:,1);
data2.data(:,3) = data2.data(:,3) ./ data2.data(:,1);
data2.data(:,1) = [];

concentrations = unique(crr(:,1));
concentrations(1) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
chesi_gitler_2012.orfs = t;
chesi_gitler_2012.data = t2;
chesi_gitler_2012.ph = strcat(phenotypes, '; ', treatments, '; ', num2str(concentrations), ' mM');

save('./chesi_gitler_2012.mat','chesi_gitler_2012');
return;

% Save data into database
dt = chesi_gitler_2012;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./chesi_gitler_2012.txt','w');
write_matrix_file(fid, chesi_gitler_2012.orfs, chesi_gitler_2012.ph, chesi_gitler_2012.data);
fclose(fid);

end

