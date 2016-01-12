%% van Voorst~Bradt, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

van_voorst_bradt_2006.pmid = 16598687;

phenotypes = {'growth [streaks on agar]'};
treatments = {'ethanol [6%]'};

% Load data
fid = fopen('./raw_data/hits_genenames.txt');
hits = textscan(fid,'%s');
hits_genenames = hits{1};
fclose(fid);

hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

van_voorst_bradt_2006.orfs = hits_orfs;
van_voorst_bradt_2006.data = hits_scores;

van_voorst_bradt_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('./van_voorst_bradt_2006.mat','van_voorst_bradt_2006');
return;

% Save data into database
dt = van_voorst_bradt_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

