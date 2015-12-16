%% Mira~Sa-Correia, 2010
function FILENAMES = code()
FILENAMES = {};
mira_sa_correia_2010.pmid = 20973990;

phenotypes = {'growth'};
treatments = {'acetic acid'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','./raw_data/List of strains tested.xlsx');
tested_orfs = tested.raw(2:end,1);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/1475-2859-9-79-s1.xlsx');
hits_genenames = data.raw(9:end,1);
hits_data = data.raw(9:end,3);

inds = find(cellfun(@isempty, hits_genenames) | cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_data(inds,:) = [];

hits_genenames = strtrim(hits_genenames);

hits_orfs = genename2orf(hits_genenames);
hits_orfs(strcmp('ace1', hits_orfs)) = {'YGL166W'};
hits_orfs(strcmp('api2', hits_orfs)) = {'YDR525W'};
hits_orfs(strcmp('brp1', hits_orfs)) = {'YGL007W'};
hits_orfs(strcmp('bud19', hits_orfs)) = {'YJL188C'};
hits_orfs(strcmp('bud26', hits_orfs)) = {'YDR241W'};
hits_orfs(strcmp('cos16', hits_orfs)) = {'YCR044C'};
hits_orfs(strcmp('cup5', hits_orfs)) = {'YEL027W'};
hits_orfs(strcmp('kem1', hits_orfs)) = {'YGL173C'};
hits_orfs(strcmp('opi8', hits_orfs)) = {'YKR035C'};
hits_orfs(strcmp('opi9', hits_orfs)) = {'YLR338W'};
hits_orfs(strcmp('rhr2', hits_orfs)) = {'YIL053W'};
hits_orfs(strcmp('see1', hits_orfs)) = {'YIL064W'};
hits_orfs(strcmp('sur4', hits_orfs)) = {'YLR372W'};
hits_orfs(strcmp('tfp1', hits_orfs)) = {'YDL185W'};
hits_orfs(strcmp('vps61', hits_orfs)) = {'YDR136C'};
hits_orfs(strcmp('vps66', hits_orfs)) = {'YPR139C'};

inds = find(strcmp('ppa1', hits_orfs));  % ambiguous genename
hits_orfs(inds) = [];
hits_data(inds,:) = [];


hits_orfs = upper(hits_orfs);
inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

hits_data(strcmp('++', hits_data)) = {-2};
hits_data(strcmp('+', hits_data)) = {-1};

hits_data = cell2mat(hits_data);
hits_data(isnan(hits_data)) = 0;

[missing, ix] = setdiff(hits_orfs, tested_orfs);    % 8 ORFs deleted
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

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'mira_sa_correia_2010.mat'],'mira_sa_correia_2010');
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

