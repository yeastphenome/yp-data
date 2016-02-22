%% van Voorst~Bradt, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

van_voorst_bradt_2006.pmid = 16598687;

phenotypes = {'growth [streaks on agar]'};
treatments = {'ethanol [6%]'};

% Load data
[FILENAMES, hits] = read_data('fopen', './raw_data/hits_genenames.txt', '%s');
hits_genenames = hits{1};

hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

van_voorst_bradt_2006.orfs = hits_orfs;
van_voorst_bradt_2006.data = hits_scores;

van_voorst_bradt_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('./van_voorst_bradt_2006.mat','van_voorst_bradt_2006');

fid = fopen('./van_voorst_bradt_2006.txt','w');
write_matrix_file(fid, van_voorst_bradt_2006.orfs, van_voorst_bradt_2006.ph, van_voorst_bradt_2006.data);
fclose(fid);

end
