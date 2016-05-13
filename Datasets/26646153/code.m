%% Choi~Basrai, 2015
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
choi_basrai_2015.pmid = 26646153;

phenotypes = {'growth'};
treatment = {'nicotinamide, 120 mM'};

%% Hit Strains

% Load file from first test
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P1');
all_data = data;
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P2');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P3');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P4');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P5');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P6');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P7');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P8');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P9');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P10');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P11');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P12');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P13');
all_data = [all_data; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS2.xlsx', 'S1 P14');
all_data = [all_data; data];

% Get the list of ORFs
strains = all_data(:,1);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(strains));
all_data(inds, :) = [];
strains(inds, :) = [];

% Clean up ORFs
strains = clean_orf(strains);

% Get data from hits
hit_data = all_data(:,6);
indx = ~cellfun(@isnumeric, hit_data);
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

% Load file from second test
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P1');
all_data2 = data;
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P2');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P3');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P4');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P5');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P6');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P7');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P8');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P9');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P10');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P11');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P12');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P13');
all_data2 = [all_data2; data];
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/supp_g3.115.022244_TableS3.xlsx', 'S2 P14');
all_data2 = [all_data2; data];

% Get the list of ORFs
strains2 = all_data2(:,1);

% Find anything that doesn't look like an ORF and remove it
inds = find(~is_orf(strains2));
all_data2(inds, :) = [];
strains2(inds, :) = [];

% Clean up ORFs
strains2 = clean_orf(strains2);

% Get data from hits
hit_data2 = all_data2(:,6);
indx = ~cellfun(@isnumeric, hit_data2);
hit_data2(indx) = {NaN};
hit_data2 = cell2mat(hit_data2);

% Combine the data and find the average
final_data = [hit_data hit_data2];
final_data = sum(final_data, 2);
final_data = 1 ./ final_data;
final_strains = strains;

% Average any repeated value
[final_strains, final_data] = grpstats(final_data, final_strains, {'gname','mean'});

% Prepare final dataset
choi_basrai_2015.orfs = final_strains;
choi_basrai_2015.ph = strcat(phenotypes, '; ', treatment);
choi_basrai_2015.data = final_data;

%% Save

save('./choi_basrai_2015.mat','choi_basrai_2015');

fid = fopen('./choi_basrai_2015.txt','w');
write_matrix_file(fid, choi_basrai_2015.orfs, choi_basrai_2015.ph, choi_basrai_2015.data);
fclose(fid);

end
