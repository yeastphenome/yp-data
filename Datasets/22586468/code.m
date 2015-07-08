%% Peyroche~Plateau, 2012
% DATA = peyroche_plateau_2012
function FILENAMES = code()
FILENAMES = {};


peyroche_plateau_2012.source = {'http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0036343.s004'};
peyroche_plateau_2012.downloaddate = {'2013-03-08'};
peyroche_plateau_2012.pmid = 22586468;
peyroche_plateau_2012.desc = {'Relative fitness defect: rf = log2(wt(Se)/mut(Se) - wt/mut + 1), where wt and mut are the generation times of the WT and mutant strains with and without Na2Se.'};

phenotypes = {'Growth, log2 ratio'};
treatments = {'Sodium selenide, 1 uM 16 h'; 'Sodium selenide, 2 uM 16 h'; 'Sodium selenide, 1 uM 27 h'; 'Sodium selenide, 2 uM 27 h'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/journal.pone.0036343.s004.xlsx', 'data');

% Get indices of the data columns
ind_data = 5:8;

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

data2.orfs = upper(data.raw(:,1));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
peyroche_plateau_2012.orfs = t;
peyroche_plateau_2012.data = t2;
peyroche_plateau_2012.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'peyroche_plateau_2012.mat'],'peyroche_plateau_2012');
return;

% Save data into database
dt = peyroche_plateau_2012;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).reporter;
    datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));


end

