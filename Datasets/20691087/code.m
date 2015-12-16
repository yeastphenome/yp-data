%% Alamgir~Golshani, 2010
function FILENAMES = code()
FILENAMES = {};

alamgir_golshani_2010.source = {'http://www.biomedcentral.com/content/supplementary/1472-6769-10-6-s1.xls'};
alamgir_golshani_2010.downloaddate = {'2013-03-06'};
alamgir_golshani_2010.pmid = 20691087;

phenotypes = {'Growth, colony size'};
treatments = {'3-AT, 22 mg/ml';'cycloheximide, 45 ng/ml';'streptomycin, 40 mg/ml';'neomycin, 5.5 mg/ml'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/1472-6769-10-6-s1.xlsx', 'Raw genome-wide data');

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

inds = setdiff(1:length(data.raw(:,1)),strmatch('Y', data.raw(:,1)));
data.raw(inds,:) = [];

data2.orfs = upper(data.raw(:, 1));
data2.data = data.raw(:, 3:14);

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Eliminate zeros
data2.data(data2.data == 0) = NaN;

% Average the 3 replicates
data3.orfs = data2.orfs;
data3.data(:,1) = nanmean(data2.data(:,1:3),2);
data3.data(:,2) = nanmean(data2.data(:,4:6),2);
data3.data(:,3) = nanmean(data2.data(:,7:9),2);
data3.data(:,4) = nanmean(data2.data(:,10:12),2);

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
alamgir_golshani_2010.orfs = t;
alamgir_golshani_2010.data = t2;
alamgir_golshani_2010.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'alamgir_golshani_2010.mat'],'alamgir_golshani_2010');
return;

% Save data into database
dt = alamgir_golshani_2010;

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

