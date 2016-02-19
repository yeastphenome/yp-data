%% Parsons~Boone, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
parsons_boone_2006.pmid = 16901791;

phenotypes = {'sensitivity'};

%% Hit Strains

% Load hit strains and tested strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/mmc8.xlsx', 'Supplementary Table 7');

% Get list of ORFs from both sets of ORFs
hit_orfs = data(4:end,1);

% Get Treatments
treatments = data(3, 2:end);

% Get rid of empty rows
inds = find(cellfun(@isnumeric, hit_orfs));
hit_orfs(inds) = [];

% Clean up ORFs
hit_orfs = clean_orf(hit_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds)); 

% Make sure all are unique
hit_orfs = unique(hit_orfs);

% Retrieve data
hit_data = data(4:end, 2:end);

% Make a matrix of NaNs
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {nan};
hit_data = cell2mat(hit_data);

% Data Calculation.
hit_data = hit_data ./ -1;

% Average any repeated value
[hit_orfs, hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% Prepare final dataset
parsons_boone_2006.orfs = hit_orfs;
parsons_boone_2006.ph = strcat(phenotypes, '; ', treatments);
parsons_boone_2006.data = hit_data;

%% Save

save('./parsons_boone_2006.mat','parsons_boone_2006');

fid = fopen('./parsons_boone_2006.txt','w');
write_matrix_file(fid, parsons_boone_2006.orfs, parsons_boone_2006.ph, parsons_boone_2006.data);
fclose(fid);

end
