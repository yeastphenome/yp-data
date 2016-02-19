%% Askree~McEachern, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};


askree_mceachern_2004.source = {'http://www.pnas.org/content/suppl/2004/05/20/0401263101.DC1/01263Table3.xls; Michael McEachern'};
askree_mceachern_2004.downloaddate = {'2013-07-17'};
askree_mceachern_2004.pmid = 15161972;

phenotypes = {'Telomere length'};
treatments = {''};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/01263Table3.xlsx', 'Sheet1');

% Get indices of the data columns
ind_data = 3:4;


% Eliminate anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric, data.raw(:,1)));
data.raw(inds,:) = [];

% Set data for inconsistent mutants to NaN
inds = find(strncmp('*Y', data.raw(:,1),2));
data.raw(inds,ind_data) = {NaN};

% Eliminate the asterisk before the ORF name
temp_orfs = char(data.raw(inds,1));
data.raw(inds,1) = cellstr(temp_orfs(:,2:end));

inds = find(~strncmp('Y', data.raw(:,1),1));
data.raw(inds,:) = [];

% Eliminate white spaces before/after ORF
data.raw(:,1) = cellfun(@strtrim, data.raw(:,1),'UniformOutput',0);

% Make sure all the data are numbers
t = data.raw(:,ind_data);
inds = find(~cellfun(@isnumeric, t));
t(inds) = {NaN};

data2.orfs = upper(data.raw(:,1));
data2.data(:,1) = nanmean(cell2mat(t),2);

% Load the tested genes
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/S.cerftpmata.xlsx', 'Sheet1');
tested_orfs = data.raw(:,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[t,t2] = grpstats(data2.data, data2.orfs, {'gname','mean'});
data3.orfs = t;
data3.data = t2;

% Put everything together
askree_mceachern_2004.orfs = tested_orfs;
askree_mceachern_2004.ph = phenotypes;
askree_mceachern_2004.data = zeros(length(tested_orfs),length(phenotypes));
[t,ind1,ind2] = intersect(data3.orfs, askree_mceachern_2004.orfs);
askree_mceachern_2004.data(ind2,:) = data3.data(ind1,:);



% Eliminate the essential genes
load essential_genes_100908;
[t,ind1,ind2] = intersect(askree_mceachern_2004.orfs, essential_genes);
askree_mceachern_2004.orfs(ind1) = [];
askree_mceachern_2004.data(ind1,:) = [];

save('./askree_mceachern_2004.mat','askree_mceachern_2004');
return;

% Save data into database
dt = askree_mceachern_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./askree_mceachern_2004.txt','w');
write_matrix_file(fid, askree_mceachern_2004.orfs, askree_mceachern_2004.ph, askree_mceachern_2004.data);
fclose(fid);

end

