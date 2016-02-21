%% Huang~Paulovich, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

huang_paulovich_2013.pmid = 23382077;

phenotypes = {'growth'};
treatments = {'MMS'};

% Load plate maps
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v4 0.xls', 'DATA');
tested_orfs = tested.raw(2:end,2);
tested_orfs(cellfun(@isnumeric, tested_orfs)) = [];

tested_orfs(ismember(tested_orfs,{'YLR287-A'})) = {'YLR287C-A'};
tested_orfs = unique(upper(clean_orf(tested_orfs)));

% Load data
[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/huang_paulovich_2013_hits.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% Create dataset
huang_paulovich_2013.orfs = tested_orfs;
huang_paulovich_2013.data = zeros(length(tested_orfs), length(treatments));
[t,ind1,ind2] = intersect(huang_paulovich_2013.orfs, hits_orfs);
huang_paulovich_2013.data(ind1,:) = hits_scores(ind2,:);

huang_paulovich_2013.ph = [strcat(phenotypes, '; ', treatments)];

save('./huang_paulovich_2013.mat','huang_paulovich_2013');

fid = fopen('./huang_paulovich_2013.txt','w');
write_matrix_file(fid, huang_paulovich_2013.orfs, huang_paulovich_2013.ph, huang_paulovich_2013.data);
fclose(fid);

end
