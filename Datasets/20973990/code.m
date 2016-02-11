%% Mira~Sa-Correia, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mira_sa_correia_2010.pmid = 20973990;

phenotypes = {'growth'};
treatments = {'acetic acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
[FILENAMES{end+1}, data.raw] = readdata('xlsread','./raw_data/1475-2859-9-79-s1.xlsx');
hits_genenames = data.raw(9:end,1);
hits_data = data.raw(9:end,3);

inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_data(inds,:) = [];

hits_genenames = cleanGenename(hits_genenames);

[hits_orfs,translated] = translate(hits_genenames);
hits_orfs(~translated) = [];
hits_data(~translated,:) = [];

hits_data(strcmp('++', hits_data)) = {-2};
hits_data(strcmp('+', hits_data)) = {-1};

hits_data = cell2mat(hits_data);
hits_data(isnan(hits_data)) = 0;

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 10 ORFs deleted
hits_orfs(ix) = [];
hits_data(ix,:) = [];


[t,t2] = grpstats(hits_data, hits_orfs, {'mean','gname'});
hits_orfs = t2;
hits_data = t;

mira_sa_correia_2010.orfs = tested_orfs;
mira_sa_correia_2010.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
mira_sa_correia_2010.data(ind2) = hits_data(ind1);

mira_sa_correia_2010.ph = strcat(phenotypes, '; ', treatments);

save('./mira_sa_correia_2010.mat','mira_sa_correia_2010');
return;

% Save data into database
dt = mira_sa_correia_2010;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));


end

