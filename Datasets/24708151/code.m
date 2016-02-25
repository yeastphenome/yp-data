%% Johnson~Wu, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
schlecht_stonge_2014.pmid = 24708151;

phenotypes = {'growth'};
treatments = {'YPD, 500 mM CuSO4'; 'YPE, 500 mM CuSO4'; 'YPG, 500 mM CuSO4'; 'YPL, 500 mM CuSO4'; 'respiratory deficient pools in YPE'};

%% Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/Table.S1.xlsx', 'Table S1');

% Get the list of ORFs
hit_strains = data(2:end, 3);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get data from hits and standardize it
hit_data = data(2:end, 8:2:16);
hit_data(strcmp('NA', hit_data)) = {NaN};
hit_data = cell2mat(hit_data);

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Prepare final dataset
schlecht_stonge_2014.orfs = hit_strains;
schlecht_stonge_2014.ph = strcat(phenotypes, '; ', treatments);
schlecht_stonge_2014.data = hit_data;

%% Save

save('./schlecht_stonge_2014.mat','schlecht_stonge_2014');

fid = fopen('./schlecht_stonge_2014.txt','w');
write_matrix_file(fid, schlecht_stonge_2014.orfs, schlecht_stonge_2014.ph, schlecht_stonge_2014.data);
fclose(fid);

end
