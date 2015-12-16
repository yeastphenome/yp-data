%% Baryshnikova~Myers, 2010
function FILENAMES = code()
FILENAMES = {};
baryshnikova_myers_2010.pmid = 21076421;

phenotypes = {'growth'};
treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/Supplementary_data_1_SMF_standard_100209.xlsx');

hits_orfs = data.raw(:,1);
hits_data = data.raw(:,2);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Eliminate the data for non-deletion strains
t = regexp(hits_orfs,'_','split');
inds = find(cellfun(@length, t) > 1);
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hits_data));
hits_data = cell2mat(hits_data);

hits_orfs = upper(strtrim(hits_orfs)); %changed second hit_orfs from tested_orfs
[t,t2] = grpstats(hits_data, hits_orfs,{'mean','gname'});

baryshnikova_myers_2010.orfs = t2;
baryshnikova_myers_2010.data = t;
baryshnikova_myers_2010.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'baryshnikova_myers_2010.mat'],'baryshnikova_myers_2010');
return;

% Save data into database
dt = baryshnikova_myers_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

