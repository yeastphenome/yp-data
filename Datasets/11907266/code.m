%% Dimmer~Westermann, 2002
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
dimmer_westermann_2002.pmid = 11907266;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(dimmer_westermann_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

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

hit_strains = tested_strains;
hit_data = zeros(length(tested_strains),1);

[~,ind1,ind2] = intersect(hit_strains, pet_mutants2);
hit_data(ind1) = -0.5;

[~,ind1,ind2] = intersect(hit_strains, pet_mutants);
hit_data(ind1) = -1;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [470];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
dimmer_westermann_2002.orfs = hit_strains;
dimmer_westermann_2002.ph = hit_data_names;
dimmer_westermann_2002.data = hit_data;
dimmer_westermann_2002.dataset_ids = hit_data_ids;

%% Save

save('./dimmer_westermann_2002.mat','dimmer_westermann_2002');

%% Print out

fid = fopen('./dimmer_westermann_2002.txt','w');
write_matrix_file(fid, dimmer_westermann_2002.orfs, dimmer_westermann_2002.ph, dimmer_westermann_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(dimmer_westermann_2002)
end

end
