%% Dos Santos~Sa-Correia, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

dos_santos_sa_correia_2011.pmid = 21960436;

phenotypes = {'growth [spot assay]'};
treatments = {'quinine [1.5-1.7 g/L]'};

% Load tested
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/List of strains tested.xlsx', 'Tabelle2');
tested_orfs = tested.raw(2:end,1);
% slow_growers = tested.raw(2:end,2);
% inds = find(~cellfun(@isnumeric, slow_growers));
% tested_orfs(inds) = [];

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));

% Load data
[FILENAMES{end+1}, hits_genenames_HS] = read_data('textread','./raw_data/hits_genenames_hs.txt', '%s');

hits_genenames_HS = clean_genename(hits_genenames_HS);
hits_orfs_HS = translate(hits_genenames_HS);
hits_scores_HS = zeros(length(hits_orfs_HS),1)-2;

[FILENAMES{end+1}, hits_genenames_S] = read_data('textread','./raw_data/hits_genenames_s.txt', '%s');

hits_genenames_S = clean_genename(hits_genenames_S);
hits_orfs_S = translate(hits_genenames_S);
hits_scores_S = zeros(length(hits_orfs_S),1)-1;

[FILENAMES{end+1}, hits_genenames_R] = read_data('textread','./raw_data/hits_genenames_r.txt', '%s');

hits_genenames_R = upper(clean_genename(hits_genenames_R));
[hits_orfs_R, translated] = translate(hits_genenames_R);
hits_orfs_R(~translated) = [];
hits_scores_R = zeros(length(hits_orfs_R),1)+1;

% Adjust the overlapping strains as follows: eliminate HS from S; eliminate
% HS ^ R and S ^ R from both
[~,~,ind2] = intersect(hits_orfs_HS, hits_orfs_S);
hits_orfs_S(ind2) = []; hits_scores_S(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_HS, hits_orfs_R);
hits_orfs_HS(ind1) = []; hits_scores_HS(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

[~,ind1,ind2] = intersect(hits_orfs_S, hits_orfs_R);
hits_orfs_S(ind1) = []; hits_scores_S(ind1) = [];
hits_orfs_R(ind2) = []; hits_scores_R(ind2) = [];

hits_orfs = [hits_orfs_HS; hits_orfs_S; hits_orfs_R];
hits_scores = [hits_scores_HS; hits_scores_S; hits_scores_R];

[missing, ix] = setdiff(hits_orfs, tested_orfs);

% Adjustments
hits_orfs(strcmp('YHR039C-A', hits_orfs)) = {'YHR039C-B'};
hits_orfs(strcmp('YBR089C-A', hits_orfs)) = {'YBR090C-A'};

dos_santos_sa_correia_2011.orfs = tested_orfs;
dos_santos_sa_correia_2011.data = zeros(length(tested_orfs), length(phenotypes));
[t,ind1,ind2] = intersect(hits_orfs, tested_orfs);
dos_santos_sa_correia_2011.data(ind2,1) = hits_scores(ind1);

dos_santos_sa_correia_2011.ph = [strcat(phenotypes, '; ', treatments)];

save('./dos_santos_sa_correia_2011.mat','dos_santos_sa_correia_2011');

fid = fopen('./dos_santos_sa_correia_2011.txt','w');
write_matrix_file(fid, dos_santos_sa_correia_2011.orfs, dos_santos_sa_correia_2011.ph, dos_santos_sa_correia_2011.data);
fclose(fid);

end
