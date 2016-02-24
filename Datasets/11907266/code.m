%% Dimmer~Westermann, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dimmer_westermann_2002.pmid = 11907266;

phenotypes = {'growth'};
treatments = {'Gly'};

% Load data
[FILENAMES, pet_mutants] = read_data('textscan', './raw_data/no_growth.txt', '%s %*[^\n]');
pet_mutants = clean_orf(pet_mutants);

inds = find(~is_orf(pet_mutants));
disp(pet_mutants(inds));
pet_mutants(inds) = [];

pet_mutants = unique(pet_mutants);

[FILENAMES, pet_mutants2] = read_data('textscan', './raw_data/poor_growth.txt', '%s %*[^\n]');
pet_mutants2 = clean_orf(pet_mutants2);

inds = find(~is_orf(pet_mutants2));
disp(pet_mutants2(inds));

pet_mutants2 = unique(pet_mutants2);

dimmer_westermann_2002.orfs = unique([pet_mutants; pet_mutants2]);
dimmer_westermann_2002.data = nan(length(dimmer_westermann_2002.orfs),1);

[~,ind1,ind2] = intersect(dimmer_westermann_2002.orfs, pet_mutants2);
dimmer_westermann_2002.data(ind1) = -0.5;

[~,ind1,ind2] = intersect(dimmer_westermann_2002.orfs, pet_mutants);
dimmer_westermann_2002.data(ind1) = -1;

dimmer_westermann_2002.ph = strcat(phenotypes, '; ', treatments);

save('./dimmer_westermann_2002.mat','dimmer_westermann_2002');

fid = fopen('./dimmer_westermann_2002.txt','w');
write_matrix_file(fid, dimmer_westermann_2002.orfs, dimmer_westermann_2002.ph, dimmer_westermann_2002.data);
fclose(fid);

end
