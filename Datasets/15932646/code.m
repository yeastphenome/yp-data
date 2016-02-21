%% Stepchenkova~Pavlov, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

stepchenkova_pavlov_2005.source = {'main PDF'};
stepchenkova_pavlov_2005.downloaddate = {'2014-03-10'};
stepchenkova_pavlov_2005.pmid = 15932646;

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/stepchenkova_pavlov_2005_hits.xlsx', 'Sheet1');

hits_genenames = data.raw(:,1);
inds = find(cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
data.raw(inds,:) = [];

hits_orfs = translate(hits_genenames);

hits_scores = cell2mat(data.raw(:,2:3));

% Normalize to WT
inds = find(strcmpi('WT', hits_orfs));
hits_scores = hits_scores ./ repmat(hits_scores(inds,:), length(hits_orfs),1);
hits_scores(inds,:) = [];
hits_orfs(inds) = [];

phenotypes = {'growth';'mutagenicity'};
treatments = {'HAP';'HAP'};


% Load tested genes
[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/mat_alpha_041902.txt', '%*s %s %*s %*s %*s %*s %*s %*s', 'delimiter', '\t');

inds = find(cellfun(@isnumeric, tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(upper(clean_orf(tested_orfs)));
tested_orfs(~is_orf(tested_orfs)) = [];

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
stepchenkova_pavlov_2005.orfs = tested_orfs;
stepchenkova_pavlov_2005.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(stepchenkova_pavlov_2005.orfs, hits_orfs);
stepchenkova_pavlov_2005.data(ind1,:) = hits_scores(ind2,:);

stepchenkova_pavlov_2005.ph = [strcat(phenotypes, '; ', treatments)];


save('./stepchenkova_pavlov_2005.mat','stepchenkova_pavlov_2005');

fid = fopen('./stepchenkova_pavlov_2005.txt','w');
write_matrix_file(fid, stepchenkova_pavlov_2005.orfs, stepchenkova_pavlov_2005.ph, stepchenkova_pavlov_2005.data);
fclose(fid);

end
