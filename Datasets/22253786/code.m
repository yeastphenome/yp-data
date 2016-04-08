%% Blackman~Nislow, 2012
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
blackman_nislow_2012.pmid = 22253786;

phenotypes = {'growth'};
treatment = {'elesclomol-Cu'};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pone.0029798.s003.xlsx', 'elescomol');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = cellfun(@(x) strtok(x, ':'), strains, 'UniformOutput', false); 
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = cell2mat(data(2:end, 3)); 

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% Prepare final dataset
blackman_nislow_2012.orfs = strains;
blackman_nislow_2012.ph = strcat(phenotypes, '; ', treatment);
blackman_nislow_2012.data = hit_data;

%% Save

save('./blackman_nislow_2012.mat','blackman_nislow_2012');

fid = fopen('./blackman_nislow_2012.txt','w');
write_matrix_file(fid, blackman_nislow_2012.orfs, blackman_nislow_2012.ph, blackman_nislow_2012.data);
fclose(fid);

end
