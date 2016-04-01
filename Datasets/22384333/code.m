%% Hoon~Nislow, 2011
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
hoon_nislow_2011.pmid = 22384333;

phenotypes = {'growth'};
treatment = {'glyoxal 80mM'};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_1.3.219_TableS3.xlsx', 'homozygous deletions');

% Get the list of ORFs
strains = data(5:end, 1);

% Clean up ORFs
strains = cellfun(@(x) strtok(x, ':'), strains, 'UniformOutput', false); 
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = -cell2mat(data(5:end, 6)); 

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% Prepare final dataset
hoon_nislow_2011.orfs = strains;
hoon_nislow_2011.ph = strcat(phenotypes, '; ', treatment);
hoon_nislow_2011.data = hit_data;

%% Save

save('./hoon_nislow_2011.mat','hoon_nislow_2011');

fid = fopen('./hoon_nislow_2011.txt','w');
write_matrix_file(fid, hoon_nislow_2011.orfs, hoon_nislow_2011.ph, hoon_nislow_2011.data);
fclose(fid);

end
