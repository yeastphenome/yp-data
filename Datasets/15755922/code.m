%% Gardocki~Lopes, 2005
% DATA = gardocki_lopes_2005
function FILENAMES = code()
FILENAMES = {};
gardocki_lopes_2005.pmid = 15755922;

phenotypes = {'expression of PIS1'};
treatments = {'glucose 2%';'glycerol 3%'};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/120_GENES_AFFECTING_PIS1.xlsx');

orfs = data.raw(8:end,2);
raw_data = data.raw(8:end,[11 14]);

inds = find(cellfun(@isempty, orfs) | cellfun(@isnumeric, orfs));
orfs(inds) = [];
raw_data(inds,:) = [];

inds = find(~strncmp('Y', orfs,1));
orfs(inds) = [];
raw_data(inds,:) = [];

% Transform symbols to numbers
% Glucose: + -> -1, - -> 0
% Glycerol: - -> 0, + -> +1

inds = find(strcmp('+', raw_data(:,1)));
raw_data(inds,1) = {-1};
inds = find(strcmp('-', raw_data(:,1)));
raw_data(inds,1) = {0};
inds = find(strcmp('+', raw_data(:,2)));
raw_data(inds,2) = {1};
inds = find(strcmp('-', raw_data(:,2)));
raw_data(inds,2) = {0};
inds = find(strcmp('NT', raw_data));
raw_data(inds) = {NaN};

raw_data = cell2mat(raw_data);

orfs = upper(strtrim(orfs));
[t,t2] = grpstats(raw_data, orfs,{'mean','gname'});

gardocki_lopes_2005.orfs = t2;
gardocki_lopes_2005.data = t;
gardocki_lopes_2005.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'gardocki_lopes_2005.mat'],'gardocki_lopes_2005');
return;

% Save data into database
dt = gardocki_lopes_2005;

% datasets = get_datasets_for_paper(dt);
%
% [~,database_ix] = sortrows(datasets.names);
% [~,ph_ix] = sort(dt.ph);
ph_ix = 1:length(dt.ph);
%
% % Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
% datasets.names(database_ix,:)
% dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, [575 576]);

end

