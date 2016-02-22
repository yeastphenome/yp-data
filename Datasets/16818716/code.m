%% Lam~Conibear, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
lam_conibear_2006.pmid = 16818716;

phenotypes = {'polytopic membrane protein trafficking'};
treatments = {''};

% Load data
[FILENAMES, C] = read_data('fopen', './raw_data/hits_genes_data.txt', '%s\t%.3f\n');

hits_genes = C{1};
hits_data = C{2};

hits_orfs = translate(hits_genes);

lam_conibear_2006.orfs = hits_orfs;
lam_conibear_2006.data = hits_data;

lam_conibear_2006.ph = strcat(phenotypes, '; ', treatments);

save('./lam_conibear_2006.mat','lam_conibear_2006');

fid = fopen('./lam_conibear_2006.txt','w');
write_matrix_file(fid, lam_conibear_2006.orfs, lam_conibear_2006.ph, lam_conibear_2006.data);
fclose(fid);

end
