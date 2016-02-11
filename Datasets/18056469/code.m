%% Xia~Flores-Rozas, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

xia_flores_rozas_2007.pmid = 18056469;

phenotypes = {'growth (colony size)'};
treatments = {'doxorubicin [20 umol/L]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = readdata('xlsread','./raw_data/Mat_a.xlsx', 'mat_a_041902');
tested_orfs = tested.raw(3:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs(strcmp('YLR287-A', tested_orfs)) = {'YLR287C-A'};
tested_orfs = unique(upper(cleanOrf(tested_orfs)));

% Load data
[FILENAMES{end+1}, DATA] = readdata('textread','./raw_data/hits_genenames.txt', '%s %s', 'delimiter', '\t');

hits_genenames = DATA{1};
hits_scores_txt = DATA{2};

hits_genenames = cleanGenename(hits_genenames);
hits_orfs = translate(hits_genenames);
hits_scores = -cellfun(@length, hits_scores_txt);

[missing, ix] = setdiff(hits_orfs, tested_orfs);


xia_flores_rozas_2007.orfs = tested_orfs;
xia_flores_rozas_2007.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
xia_flores_rozas_2007.data(ind2,1) = hits_scores(ind1);

xia_flores_rozas_2007.ph = [strcat(phenotypes, '; ', treatments)];

save('./xia_flores_rozas_2007.mat','xia_flores_rozas_2007');
return;

% Save data into database
dt = xia_flores_rozas_2007;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

