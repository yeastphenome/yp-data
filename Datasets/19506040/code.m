%% Burston~Conibear, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
burston_conibear_2009.pmid = 19506040;

phenotypes = {'endocytosis (MatA)';'endocytosis (MatAlpha)'};
treatments = {''};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/JCB_200811116_TS1.xlsx', 'TableS1');

hits_orfs = data.raw(6:end,2);
hits_data_a = data.raw(6:end,4);
hits_data_alpha = data.raw(6:end,5);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data_a(inds) = [];
hits_data_alpha(inds) = [];

burston_conibear_2009.orfs = hits_orfs;
burston_conibear_2009.data = cell2mat([hits_data_a hits_data_alpha]);

burston_conibear_2009.ph = strcat(phenotypes, '; ', treatments);

save('./burston_conibear_2009.mat','burston_conibear_2009');

fid = fopen('./burston_conibear_2009.txt','w');
write_matrix_file(fid, burston_conibear_2009.orfs, burston_conibear_2009.ph, burston_conibear_2009.data);
fclose(fid);

end
