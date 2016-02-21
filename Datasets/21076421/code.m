%% Baryshnikova~Myers, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
baryshnikova_myers_2010.pmid = 21076421;

phenotypes = {'growth'};
treatments = {'standard'};

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/Supplementary_data_1_SMF_standard_100209.xlsx');

hits_orfs = data.raw(:,1);
hits_data = data.raw(:,2);

inds = find(cellfun(@isempty, hits_orfs) | cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

% Eliminate the data for non-deletion strains
t = regexp(hits_orfs,'_','split');
inds = find(cellfun(@length, t) > 1);
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds,:) = [];

inds = find(~cellfun(@isnumeric, hits_data));
hits_data = cell2mat(hits_data);

hits_orfs = upper(strtrim(hits_orfs)); %changed second hit_orfs from tested_orfs
[t,t2] = grpstats(hits_data, hits_orfs,{'mean','gname'});

baryshnikova_myers_2010.orfs = t2;
baryshnikova_myers_2010.data = t;
baryshnikova_myers_2010.ph = strcat(phenotypes, '; ', treatments);

save('./baryshnikova_myers_2010.mat','baryshnikova_myers_2010');

fid = fopen('./baryshnikova_myers_2010.txt','w');
write_matrix_file(fid, baryshnikova_myers_2010.orfs, baryshnikova_myers_2010.ph, baryshnikova_myers_2010.data);
fclose(fid);

end
