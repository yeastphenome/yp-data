%% Yoshikawa~Shimizu, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

yoshikawa_shimizu_2009.pmid = 19054128;
yoshikawa_shimizu_2009.desc = {'The final values are treated divided by untreated.'};

phenotypes = {'Growth, exponential growth rate'};
treatments = {'UNT';'EtOH, 5%';'EtOH, 8%';'NaCl, 1 M'};

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/FYR_456_sm_tableS1.xlsx', 'data');

data2.orfs = upper(data.raw(3:end, 1));

ind_data = [3 4 6 7 9 10 12 13];
data2.data = data.raw(3:end, ind_data);

% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

inds = setdiff(1:length(data2.orfs),strmatch('Y', data2.orfs));
data2.orfs(inds) = [];
data2.data(inds,:) = [];

% Make sure all the data are numbers
inds = find(cellfun(@isnumeric, data2.data)==0);
data2.data(inds) = {NaN};
data2.data = cell2mat(data2.data);

% Average the 2 replicates
data3.orfs = data2.orfs;
data3.data(:,1) = nanmean(data2.data(:,1:2),2);
data3.data(:,2) = nanmean(data2.data(:,3:4),2);
data3.data(:,3) = nanmean(data2.data(:,5:6),2);
data3.data(:,4) = nanmean(data2.data(:,7:8),2);

% Normalize by UNTREATED sample
data3.data(:,2) = data3.data(:,2)./data3.data(:,1);
data3.data(:,3) = data3.data(:,3)./data3.data(:,1);
data3.data(:,4) = data3.data(:,4)./data3.data(:,1);
data3.data(:,1) = [];
treatments(1) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data3.data, data3.orfs, {'gname','mean'});
yoshikawa_shimizu_2009.orfs = t;
yoshikawa_shimizu_2009.data = t2;
yoshikawa_shimizu_2009.ph = strcat(phenotypes, {'; '}, treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'yoshikawa_shimizu_2009.mat'],'yoshikawa_shimizu_2009');
return;

% Save data into database
dt = yoshikawa_shimizu_2009;

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
datasets_names(database_ix([2 3 1]),:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets_ids(database_ix([2 3 1])));

end

