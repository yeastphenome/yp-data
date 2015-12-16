%% Hartman~Tippery, 2004
function FILENAMES = code()
FILENAMES = {};

hartman_tippery_2004.pmid = 15239834;
hartman_tippery_2004.desc = {'The loaded values are Z-scores with respect to the wild-type in a given condition (UNT or HU). The final values (normalized to UNT) are the difference between treated and untreated z-scores.'};

phenotypes = {'Growth, AUGC'};
treatments = {'UNT';'HU, 50 mM';'HU, 150 mM'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/gb-2004-5-7-r49-s7.xlsx', 'data');

ind_orf = strmatch('ORF', data.raw(1,:));
data2.orfs = data.raw(2:end, ind_orf);

ind_data1 = strmatch('No HU- Growth Index', data.raw(1,:));
ind_data2 = strmatch('50mM HU Growth Index', data.raw(1,:));
ind_data3 = strmatch('150 mM HU Growth Index', data.raw(1,:));

data2.data = data.raw(2:end, [ind_data1 ind_data2 ind_data3]);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

data2.orfs = upper(data2.orfs);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Normalize by UNTREATED sample
data2.data(:,2) = data2.data(:,2) - data2.data(:,1);
data2.data(:,3) = data2.data(:,3) - data2.data(:,1);
data2.data(:,1) = [];
treatments(1) = [];


% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
hartman_tippery_2004.orfs = t;
hartman_tippery_2004.data = t2;
hartman_tippery_2004.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'hartman_tippery_2004.mat'],'hartman_tippery_2004');
return;

% Save data into database
dt = hartman_tippery_2004;

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

