%% Reiner~Schneiter, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

reiner_schneiter_2006.pmid = 16251356;

phenotypes = {'growth'};
treatments = {'hypoxia'};

% Load tested genes
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/BY4741-MATa COLLECTION.xls', 'chr11_1yes');

tested_orfs = data.raw(2:end,2);

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

% Eliminate anything that doesn't look like an ORF
inds = find(~strncmp('Y', tested_orfs,1));
tested_orfs(inds) = [];

tested_orfs = unique(upper(tested_orfs));

% Load hits
[FILENAMES{end+1}, hits] = read_data('textread','./raw_data/reiner_scheiter_2006_hits.txt', '%s');

% This list of ORFs is lacking the last character (published that way), so
% we have to match it to the list of tested strains.
for i = 1 : length(hits)
inds = find(strncmp(hits{i}, tested_orfs, length(hits{i})));
if length(inds) == 1
hits_orfs(i) = tested_orfs(inds);
else
fprintf('%s\t%d\n', hits{i}, length(inds));
end
end

% Two ORFs (YBR039W and YNL243W) could not be found in the list of tested
% strains. So we have to remove them.
hits_orfs(strncmp('YBR039', hits, length('YBR039')) | strncmp('YNL243', hits, length('YNL243'))) = [];

% Eliminate white spaces before/after ORF
hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

% Assign a -1 to all the strains.
hits_scores = zeros(length(hits_orfs),1)-1;

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
reiner_schneiter_2006.orfs = tested_orfs;
reiner_schneiter_2006.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(reiner_schneiter_2006.orfs, hits_orfs);
reiner_schneiter_2006.data(ind1,:) = hits_scores(ind2,:);

reiner_schneiter_2006.ph = [strcat(phenotypes, '; ', treatments)];

save('./reiner_schneiter_2006.mat','reiner_schneiter_2006');

fid = fopen('./reiner_schneiter_2006.txt','w');
write_matrix_file(fid, reiner_schneiter_2006.orfs, reiner_schneiter_2006.ph, reiner_schneiter_2006.data);
fclose(fid);

end
