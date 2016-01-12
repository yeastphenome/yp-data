%% Chavel~Cullen, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

chavel_cullen_2010.source = {'http://www.plosgenetics.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pgen.1000883.s011'};
chavel_cullen_2010.downloaddate = {'2013-03-06'};
chavel_cullen_2010.pmid = 20333241;

phenotypes = {'Msb2p secretion'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/journal.pgen.1000883.s011.xlsx', 'Complete Screen');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,2)));
data.raw(inds,:) = [];

% A couple of manual fixes
data.raw(strcmp('YLR287-A', data.raw(:,2)),2) = {'YLR287C-A'};
data.raw(strcmp('YML009c', data.raw(:,2)),2) = {'YML009C'};

inds = find(cellfun(@isempty, regexp(data.raw(:,2), 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*')));
data.raw(inds,:) = [];

% Eliminate false positives
inds = find(~cellfun(@isnumeric, data.raw(:,7)));
data.raw(inds,6) = {NaN};

data2.orfs = upper(data.raw(:,2));
data2.data = data.raw(:,6);

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, data2.data));
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Eliminate zeros
data2.data(data2.data == 0) = NaN;

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
chavel_cullen_2010.orfs = t;
chavel_cullen_2010.data = t2;
chavel_cullen_2010.ph = phenotypes;

save('./chavel_cullen_2010.mat','chavel_cullen_2010');
return;

% Save data into database
dt = chavel_cullen_2010;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
datasets_names{i,2} = datasets(i).shortname;
datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

end

