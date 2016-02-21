%% Dakshinamurthy~Garfinkel, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dakshinamurthy_garfinkel_2010.pmid = 20498295;

phenotypes = {'decreased Ty1 transposon mobility'};
treatments = {''};

load /Datasets/18202368/nyswaner_garfinkel_2008;
dakshinamurthy_garfinkel_2010.orfs = nyswaner_garfinkel_2008.orfs;      % Same screen, positive & negative phenotypes analyzed separately
dakshinamurthy_garfinkel_2010.data = zeros(length(dakshinamurthy_garfinkel_2010.orfs),1);

% Load data
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TableS1-2.xlsx', 'Sheet2');
hits_orfs = data.raw(:,2);
hits_data = data.raw(:,3);

inds = find(cellfun(@isempty, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

inds = find(cellfun(@isnumeric, hits_orfs));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_orfs = cellfun(@strtrim, hits_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', hits_orfs,1));
hits_orfs(inds) = [];
hits_data(inds) = [];

hits_data = cellfun(@length, hits_data);

missing = setdiff(hits_orfs, dakshinamurthy_garfinkel_2010.orfs);

% Adding 4 ORFs to the list of tested
dakshinamurthy_garfinkel_2010.orfs = [dakshinamurthy_garfinkel_2010.orfs; missing];

[~,ind1,ind2] = intersect(dakshinamurthy_garfinkel_2010.orfs, hits_orfs);
dakshinamurthy_garfinkel_2010.data(ind1,1) = hits_data(ind2);

dakshinamurthy_garfinkel_2010.ph = strcat(phenotypes, '; ', treatments);

save('./dakshinamurthy_garfinkel_2010.mat','dakshinamurthy_garfinkel_2010');

fid = fopen('./dakshinamurthy_garfinkel_2010.txt','w');
write_matrix_file(fid, dakshinamurthy_garfinkel_2010.orfs, dakshinamurthy_garfinkel_2010.ph, dakshinamurthy_garfinkel_2010.data);
fclose(fid);

end
