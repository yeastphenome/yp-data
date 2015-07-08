%% Warringer~Blomberg, 2003
% DATA = warringer_blomberg_2003
function FILENAMES = code()
FILENAMES = {};

warringer_blomberg_2003.pmid = 14676322;

phenotypes = {'Growth, lag phase';'Growth, exponential growth rate';'Growth, saturation level'};
treatments = {'NaCl, 0.85 M'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/ORIG130305_LPI NaCl.xlsx', 'LPI');

data2.orfs = upper(data.raw(5:end,1));
data2.data = cell2mat(data.raw(5:end,2:4));

% Eliminate anything that doesn't look like an ORF
inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
warringer_blomberg_2003.orfs = t;
warringer_blomberg_2003.data = t2;
warringer_blomberg_2003.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'warringer_blomberg_2003.mat'],'warringer_blomberg_2003');
return;

% Save data into database
datasets = get_datasets_for_paper(warringer_blomberg_2003);
datasets_ids = zeros(length(datasets),1);
datasets_names = cell(length(datasets),3);
for i = 1 : length(datasets)
    datasets_ids(i,1) = datasets(i).id;
    datasets_names{i,1} = datasets(i).name;
    datasets_names{i,2} = datasets(i).shortname;
    datasets_names{i,3} = datasets(i).condition_dose;
end

[~,database_ix] = sortrows(datasets_names,[1 2 3]);
[~,ph_ix] = sort(warringer_blomberg_2003.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.

insert_data_into_db(warringer_blomberg_2003, ph_ix, datasets_ids(database_ix));


end

