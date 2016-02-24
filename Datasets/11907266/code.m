%% Dimmer~Westermann, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dimmer_westermann_2002.pmid = 11907266;

phenotypes = {'growth'};
treatments = {'Gly'};

%% Load tested strains
[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/HOMOZYGOUS DIPLOID 1+2 ResGen.txt', '%s %s %s %*[^\n]');
tested_strains = clean_orf(data{2});

tested_strains(strcmp('YMR41W', tested_strains)) = {'YMR241W'}; % guessed based on the position of the ORF in the file
tested_strains(strcmp('YELOO1C', tested_strains)) = {'YEL001C'};

inds2 = find(strcmp('YDL', tested_strains));
tested_strains(inds2) = strcat(tested_strains(inds2), data{3}(inds2));

inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));

tested_strains(inds) = [];
tested_strains = unique(tested_strains);


%% Load data
[FILENAMES{end+1}, pet_mutants] = read_data('textscan', './raw_data/no_growth.txt', '%s %*[^\n]');
pet_mutants = clean_orf(pet_mutants);

inds = find(~is_orf(pet_mutants));
disp(pet_mutants(inds));
pet_mutants(inds) = [];

pet_mutants = unique(pet_mutants);

[missing,~] = setdiff(pet_mutants, tested_strains); % none missing.


[FILENAMES{end+1}, pet_mutants2] = read_data('textscan', './raw_data/poor_growth.txt', '%s %*[^\n]');
pet_mutants2 = clean_orf(pet_mutants2);

inds = find(~is_orf(pet_mutants2));
disp(pet_mutants2(inds));

pet_mutants2 = unique(pet_mutants2);

[missing, ~] = setdiff(pet_mutants2, tested_strains);   % none missing

dimmer_westermann_2002.orfs = tested_strains;
dimmer_westermann_2002.data = zeros(length(tested_strains),1);

[~,ind1,ind2] = intersect(tested_strains, pet_mutants2);
dimmer_westermann_2002.data(ind1) = -0.5;

[~,ind1,ind2] = intersect(tested_strains, pet_mutants);
dimmer_westermann_2002.data(ind1) = -1;

dimmer_westermann_2002.ph = strcat(phenotypes, '; ', treatments);

save('./dimmer_westermann_2002.mat','dimmer_westermann_2002');

fid = fopen('./dimmer_westermann_2002.txt','w');
write_matrix_file(fid, dimmer_westermann_2002.orfs, dimmer_westermann_2002.ph, dimmer_westermann_2002.data);
fclose(fid);

end
