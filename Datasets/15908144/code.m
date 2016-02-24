%% Luban~Schmidt, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
luban_schmidt_2005.pmid = 15908144;

phenotypes = {'growth';'growth (after introduction of intronless mtDNA)'};
treatments = {'Gly'};

% Load tested
[FILENAMES, tested_orfs] = read_data('textscan', './raw_data/list_of_used_knockouts_PhD_Thesis_Luban.txt', '%s\n');
tested_orfs = clean_orf(tested_orfs);

inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));

tested_orfs(strcmp('YBL098\V', tested_orfs)) = {'YBL098W'};
tested_orfs(strcmp('YDR07SW', tested_orfs)) = {'YDR075W'};
tested_orfs(strcmp('YDR27SW', tested_orfs)) = {'YDR275W'};
tested_orfs(strcmp('YDR51SW', tested_orfs)) = {'YDR515W'};
tested_orfs(strcmp('YDRS41C', tested_orfs)) = {'YDR541C'};
tested_orfs(strcmp('YEL0I6C', tested_orfs)) = {'YEL016C'};
tested_orfs(strcmp('YIIL016C', tested_orfs)) = {'YHL016C'};
tested_orfs(strcmp('YIIL017W', tested_orfs)) = {'YHL017W'};
tested_orfs(strcmp('YHL0I9C', tested_orfs)) = {'YHL019C'};
tested_orfs(strcmp('YJR09JC', tested_orfs)) = {'YJR091C'};
tested_orfs(strcmp('YNL09SC', tested_orfs)) = {'YNL095C'};
tested_orfs(strcmp('YPLOI8W', tested_orfs)) = {'YPL018W'};
tested_orfs(strcmp('YPL07LC', tested_orfs)) = {'YPL071C'};

tested_orfs = unique(tested_orfs);

% Load data
[FILENAMES, pet_mutants] = read_data('textscan', './raw_data/list_of_pet_mutants.txt', '%s');
pet_mutants = clean_orf(pet_mutants);

inds = find(~is_orf(pet_mutants));
disp(pet_mutants(inds));

pet_mutants = unique(pet_mutants);

[FILENAMES, intron_mutants] = read_data('textscan', './raw_data/list_of_intron_def_mutants.txt', '%s');
intron_mutants = clean_orf(intron_mutants);

inds = find(~is_orf(intron_mutants));
disp(intron_mutants(inds));

intron_mutants = unique(intron_mutants);

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
