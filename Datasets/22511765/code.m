%% Kim~Cunningham, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

kim_cunningham_2012.pmid = 22511765;

phenotypes = {'Cell death, frequency of dead cells'};
treatments = {'tunicamycin, 2.5 ug/ml (Z-score ln)'; 'tunicamycin + FK506, 2.5 ug/ml + 1 ug/ml (Z-score ln)'; 'dithiothreitol, 4 mM (Z-score ln)'; 'dithiothreitol, 4 mM, + FK506, 1 ug/ml (Z-score ln)'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/jbc.M112.363390-1.xlsx', 'sup table 1');

% Get indices of the data columns
ind_data = [14:15 17:18];

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,5)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,5)),strmatch('Y', data.raw(:,5)));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,5) = cellfun(@strtrim, data.raw(:,5),'UniformOutput',0);

data2.orfs = upper(data.raw(:,5));
data2.data = data.raw(:,ind_data);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
kim_cunningham_2012.orfs = t;
kim_cunningham_2012.data = t2;
kim_cunningham_2012.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'kim_cunningham_2012.mat'],'kim_cunningham_2012');
return;

% Save data into database
dt = kim_cunningham_2012;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
if isempty(datasets(i).reporter)
datasets_names{i,2} = '';
else
datasets_names{i,2} = datasets(i).reporter;
end
datasets_names{i,3} = datasets(i).conditionset;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix([2 1 4 3]),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix(([2 1 4 3]))));

end

