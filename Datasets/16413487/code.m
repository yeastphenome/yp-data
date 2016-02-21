%% Alto~Dixon, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alto_dixon_2006.pmid = 16413487;

% Part 1
phenotypes = {'growth'};
treatments = {'IpgB2 effector protein'};

% Load data
[FILENAMES{end+1}, hits_gn] = read_data('textread','./raw_data/hits.txt', '%s');
hits_data = ones(size(hits_gn));

hits_orfs = translate(hits_gn);

% Load tested
[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/FG_array_genes.txt', '%s');

[missing, ix] = setdiff(hits_orfs, tested_orfs);

alto_dixon_2006.orfs = tested_orfs;
alto_dixon_2006.data = zeros(length(tested_orfs),1);

[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
alto_dixon_2006.data(ind2,:) = hits_data(ind1,:);

alto_dixon_2006.ph = strcat(phenotypes, '; ', treatments);

save('./alto_dixon_2006.mat','alto_dixon_2006');

fid = fopen('./alto_dixon_2006.txt','w');
write_matrix_file(fid, alto_dixon_2006.orfs, alto_dixon_2006.ph, alto_dixon_2006.data);
fclose(fid);

end
