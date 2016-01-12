%% Bleackley~MacGillivray, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


bleackley_macgillivray_2011.source = {'Mark Bleackley, Ross MacGillivray'};
bleackley_macgillivray_2011.notes = {'File received from the authors because online SOM didn''t have ORFs, only genenames.'};
bleackley_macgillivray_2011.downloaddate = {'2013-03-11'};
bleackley_macgillivray_2011.pmid = 21212869;

phenotypes = {'Growth, colony size'};
treatments = {'Iron, Fe(NH4)2(SO4)2 (10 mM)';'Copper, CuCl2 (7 mM)';'Manganese, MnCl2 (4 mM)';'Nickel, NiCl2 (3 mM)';'Zinc, ZnCl2 (7 mM)'; 'Cobalt, CoCl2 (2.5 mM)'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/metallomicsbleackley raw data.xls', 'rawdata');

% Get indices of the data columns
ind_data = 3:2:13;


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
bleackley_macgillivray_2011.orfs = t;
bleackley_macgillivray_2011.data = t2;
bleackley_macgillivray_2011.ph = strcat(phenotypes, {'; '}, treatments);

save('./bleackley_macgillivray_2011.mat','bleackley_macgillivray_2011');
return;

% Save data into database
dt = bleackley_macgillivray_2011;

datasets = get_datasets_for_paper(dt);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
datasets_ids(i,1) = datasets(i).id;
datasets_names{i,1} = datasets(i).name;
datasets_names{i,2} = datasets(i).short_name;
datasets_names{i,3} = datasets(i).dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets_names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix));

end

