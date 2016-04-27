%% Vahey~Voldman, 2013
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
vahey_voldman_2013.pmid = 23661198;

phenotypes = {'conductivity'};
treatment = {'electrical frequency of 300 kHz', 'electrical frequency of 10 MHz'};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/ESI_DATA_SET.txt', '%s %s %s %s %f %f %f %f');

% Get the list of ORFs
strains = data{1};

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = data{3};
hit_data(:,2) = data{4};
indx = ~cellfun(@isnumeric, hit_data);
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% Prepare final dataset
vahey_voldman_2013.orfs = strains;
vahey_voldman_2013.ph = strcat(phenotypes, '; ', treatment);
vahey_voldman_2013.data = hit_data;

%% Save

save('./vahey_voldman_2013.mat','vahey_voldman_2013');

fid = fopen('./vahey_voldman_2013.txt','w');
write_matrix_file(fid, vahey_voldman_2013.orfs, vahey_voldman_2013.ph, vahey_voldman_2013.data);
fclose(fid);

end
