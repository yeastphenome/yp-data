%% Sambade~Kane, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

sambade_kane_2005.pmid = 15937126;

phenotypes = {'growth (colony size)'};
treatments = {'pH [7.5] CaCl2 [60 mM]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/ResGen 384 well set 14 plates.xlsx', 'ResGen MATa -384');
tested_orfs = tested.raw(4:end,2);
inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs(ismember(tested_orfs,{'YLR287-A'})) = {'YLR287C-A'};

tested_orfs = unique(upper(clean_orf(tested_orfs)));

% Load data
[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/hits_genenames.txt', '%s');
hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmpi('YHR039C-A', hits_orfs)) = {'YHR039C-B'};

[missing, ix] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];   % 2 orfs to be added

sambade_kane_2005.orfs = tested_orfs;
sambade_kane_2005.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
sambade_kane_2005.data(ind2,1) = hits_scores(ind1);

sambade_kane_2005.ph = [strcat(phenotypes, '; ', treatments)];

save('./sambade_kane_2005.mat','sambade_kane_2005');

fid = fopen('./sambade_kane_2005.txt','w');
write_matrix_file(fid, sambade_kane_2005.orfs, sambade_kane_2005.ph, sambade_kane_2005.data);
fclose(fid);

end
