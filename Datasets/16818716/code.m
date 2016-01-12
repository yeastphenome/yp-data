%% Lam~Conibear, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
lam_conibear_2006.pmid = 16818716;

phenotypes = {'polytopic membrane protein trafficking'};
treatments = {''};

% Load data
fid = fopen('./raw_data/hits_genes_data.txt');
C = textscan(fid, '%s\t%.3f\n');
fclose(fid);

hits_genes = C{1};
hits_data = C{2};

hits_orfs = translate(hits_genes);

lam_conibear_2006.orfs = hits_orfs;
lam_conibear_2006.data = hits_data;

lam_conibear_2006.ph = strcat(phenotypes, '; ', treatments);

save('./lam_conibear_2006.mat','lam_conibear_2006');
return;

% Save data into database
dt = lam_conibear_2006;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[4 1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

end

