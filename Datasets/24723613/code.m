%% Lee~Giaever, 2014
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
lee_giaever_2014.pmid = 24723613;

phenotypes = {'growth'};
treatment = {'standard'};

%% Hit Strains

% Load text file
[FILENAMES{end+1}, data] = read_data('read_matrix_file', './raw_data/fitness_defect_matrix_hom.txt',1,1);

% Get Phenotypes
data.labels_col = cellfun(@strtrim, data.labels_col, 'Un', 0);
[FILENAMES{end+1}, compound] = read_data('xlsread','./raw_data/1250217s1.xlsx', 'compound library');

% Clean up the ORFs
data.labels_row = clean_orf(data.labels_row);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(data.labels_row));
disp(data.labels_row(inds));

% Match Genenames with IDs
[~, ind1, ind2] = intersect(data.labels_col, compound(:, 1));
data.labels_col_name(ind1) = compound(ind2, 2);

% Prepare final dataset
lee_giaever_2014.orfs = data.labels_row;
lee_giaever_2014.ph = strcat(phenotypes, '; ', data.labels_col_name);
lee_giaever_2014.data = data.data;

%% Save
save('./lee_giaever_2014.mat','lee_giaever_2014');

fid = fopen('./lee_giaever_2014.txt','w');
write_matrix_file(fid, lee_giaever_2014.orfs, lee_giaever_2014.ph, lee_giaever_2014.data);
fclose(fid);
end
