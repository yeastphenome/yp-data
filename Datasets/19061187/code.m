%% Matsufuji~Nakagawa, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

matsufuji_nakagawa_2008.pmid = 19061187;

phenotypes = {'growth (colony size)'};
treatments = {'acetaldehyde [35 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/S.c 5000.xlsx', 'remake');
tested_orfs = tested.raw(3:end,3);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];
tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];
tested_orfs = unique(upper(tested_orfs));

% Load data
[FILENAMES{end+1}, hits_orfs] = read_data('textread','./raw_data/hits_orfs.txt', '%s');

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);
hits_orfs = unique(upper(hits_orfs));

% Correction: ORF YHR041C was in the hit list twice. One instance should
% have been YDR378C.

hits_orfs = [hits_orfs; {'YDR378C'}];

hits_scores = -ones(length(hits_orfs),1);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_scores(inds) = [];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
tested_orfs = [tested_orfs; missing];   % 1 orfs to be added

matsufuji_nakagawa_2008.orfs = tested_orfs;
matsufuji_nakagawa_2008.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
matsufuji_nakagawa_2008.data(ind2,1) = hits_scores(ind1);

matsufuji_nakagawa_2008.ph = [strcat(phenotypes, '; ', treatments)];

save('./matsufuji_nakagawa_2008.mat','matsufuji_nakagawa_2008');

fid = fopen('./matsufuji_nakagawa_2008.txt','w');
write_matrix_file(fid, matsufuji_nakagawa_2008.orfs, matsufuji_nakagawa_2008.ph, matsufuji_nakagawa_2008.data);
fclose(fid);

end
