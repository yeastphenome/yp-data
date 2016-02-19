%% Johnson~Wu, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
johnson_wu_2015.pmid = 26688044;

phenotypes = {'growth'};
treatments = {'NaAsO2, 0.2 mM'; 'NaAsO2, 0.4 mM'};

%% Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/c5mt00261c1.xlsx', 'Sheet1');

% Get the list of ORFs
hit_strains = data(4:end, 1);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get data from hits
fid = fopen('./johnson_wu_2015.txt','w');
write_matrix_file(fid, johnson_wu_2015.orfs, johnson_wu_2015.ph, johnson_wu_2015.data);
fclose(fid);

hit_data = cell2mat(data(4:end, 2:7)); 
wt_data = cell2mat(data(3,2:7));

% Standardize the data based on 2 controls: water and wildtype
water = hit_data(:,2)-hit_data(:,1);
wt_water = wt_data(2)-wt_data(1);
water = water ./ wt_water;

concentration_one = hit_data(:,4)-hit_data(:,3);
wt_concentration_one = wt_data(4)-wt_data(3);
concentration_one = concentration_one ./ wt_concentration_one;
concentration_one = concentration_one ./ water;

concentration_two = hit_data(:,6)-hit_data(:,5);
wt_concentration_two = wt_data(6)-wt_data(5);
concentration_two = concentration_two ./ wt_concentration_two;
concentration_two = concentration_two ./ water;

% Add to make final data
hit_data = [concentration_one concentration_two];

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Prepare final dataset
johnson_wu_2015.orfs = hit_strains;
johnson_wu_2015.ph = strcat(phenotypes, '; ', treatments);
johnson_wu_2015.data = hit_data;

%% Save

save('./johnson_wu_2015.mat','johnson_wu_2015');

