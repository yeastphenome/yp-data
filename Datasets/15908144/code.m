%% Luban~Schmidt, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
luban_schmidt_2005.pmid = 15908144;

phenotypes = {'petite';'mtDNA intron slicing'};
treatments = {'Gly'};

% Load tested
fid = fopen('./raw_data/list_of_used_knockouts_PhD_Thesis_Luban.txt');
C = textscan(fid, '%s\n');
fclose(fid);
tested_orfs = unique(upper(C{1}));

% Load data
fid = fopen('./raw_data/list_of_pet_mutants.txt');
C = textscan(fid,'%s\n');
fclose(fid);
pet_mutants = unique(upper(C{1}));

fid = fopen('./raw_data/list_of_intron_def_mutants.txt');
C = textscan(fid,'%s\n');
fclose(fid);
intron_mutants = unique(upper(C{1}));

inds = find(cellfun(@isempty, tested_orfs));

inds = find(cellfun(@isnumeric, tested_orfs));


tested_orfs = cellfun(@strtrim, tested_orfs,'UniformOutput',0);

inds = find(~strncmp('Y', tested_orfs,1));

luban_schmidt_2005.orfs = tested_orfs;
luban_schmidt_2005.data = nan(length(tested_orfs),2);

missing = setdiff(pet_mutants, tested_orfs);

[~,ind1,ind2] = intersect(tested_orfs, pet_mutants);
luban_schmidt_2005.data(ind1,1) = 1;
luban_schmidt_2005.data(isnan(luban_schmidt_2005.data(:,1)),1) = 0;

luban_schmidt_2005.data(ind1,2) = 0;
[~,ind1,ind2] = intersect(tested_orfs, intron_mutants);
luban_schmidt_2005.data(ind1,2) = -1;

luban_schmidt_2005.ph = strcat(phenotypes, '; ', treatments);

save('./luban_schmidt_2005.mat','luban_schmidt_2005');

fid = fopen('./luban_schmidt_2005.txt','w');
write_matrix_file(fid, luban_schmidt_2005.orfs, luban_schmidt_2005.ph, luban_schmidt_2005.data);
fclose(fid);

end
