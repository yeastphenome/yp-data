%% Aouida~Ramotar, 2004
% DATA = aouida_ramotar_2004
function FILENAMES = code()
FILENAMES = {};
aouida_ramotar_2004.pmid = 14871844;

phenotypes = {'growth'};
treatments = {'bleomycin'};

% Load tested
[FILENAMES{end+1}, tested.raw] = dataread('xlsread','raw_data/HU haploid.xlsx');
tested_orfs = tested.raw(4:end,5);

inds = find(cellfun(@isempty, tested_orfs) | cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(strtrim(upper(tested_orfs)));

inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

% Load data
[FILENAMES{end+1}, data.raw] = dataread('xlsread','raw_data/CAN_2-1-04_Aouida.xlsx');
hits_orfs = data.raw(:,1);
hits_orfs = strtrim(upper(hits_orfs));



hits_data = strtrim(data.raw(:,2));

% Data conversion
hits_data(strcmp('R+++', hits_data)) = {9};     % (>500-fold more resistant -> log2(500) ~ 9)
hits_data(strcmp('R++++', hits_data)) = {10};
hits_data(strcmp('S+', hits_data)) = {-1};
hits_data(strcmp('S++', hits_data)) = {-2};
hits_data(strcmp('S+++', hits_data)) = {-3};
hits_data(strcmp('S++++', hits_data)) = {-4};

hits_data = cell2mat(hits_data);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];


aouida_ramotar_2004.orfs = tested_orfs;
aouida_ramotar_2004.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
aouida_ramotar_2004.data(ind2) = hits_data(ind1);

aouida_ramotar_2004.ph = strcat(phenotypes, '; ', treatments);

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'aouida_ramotar_2004.mat'],'aouida_ramotar_2004');
return;

% Save data into database
dt = aouida_ramotar_2004;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

