%% Hancock~Lopes, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hancock_lopes_2006.pmid = 16582425;

phenotypes = {'Opi-'};
treatments = {''};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/mat_alpha_061101.xlsx', 'mat_alpha_061101');
tested_orfs = tested.raw(4:end,2);

inds = cellfun(@isnumeric, tested_orfs);
tested_orfs(inds) = [];

tested_orfs = clean_orf(tested_orfs);
tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};

tested_orfs = unique(tested_orfs);

% Load data
fid = fopen('./raw_data/hits_genenames.txt');
hits_genenames = textscan(fid, '%s');
hits_genenames = hits_genenames{1};
fclose(fid);

[hits_orfs, translated] = translate(hits_genenames);
hits_orfs(~translated) = [];

hits_orfs = unique(hits_orfs);

hits_scores = ones(length(hits_orfs),1);

hancock_lopes_2006.orfs = tested_orfs;
hancock_lopes_2006.data = zeros(length(tested_orfs), length(phenotypes));
[~,ind1,ind2] = intersect(hits_orfs, tested_orfs);
hancock_lopes_2006.data(ind2,1) = hits_scores(ind1,1);

hancock_lopes_2006.ph = strcat(phenotypes, '; ', treatments);

save('./hancock_lopes_2006.mat','hancock_lopes_2006');

fid = fopen('./hancock_lopes_2006.txt','w');
write_matrix_file(fid, hancock_lopes_2006.orfs, hancock_lopes_2006.ph, hancock_lopes_2006.data);
fclose(fid);

end
