%% Lis~Bobek, 2013
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
lis_bobek_2013.pmid = 23208710;

phenotypes = {'sensitivity'};
treatments = {'MUC7 12-mer, 10 ?M'; 'histatin 12-mer, 20 ?M'; 'cathelicidin KR20, 10 ?M'; 'peptide w/ lactoferricin amino acids 1 to 11, 12 ?M'};

%% Hit Strains

% Load hit strains
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/AAC.01439-12_zac999101582so1.xlsx', 'all data');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = cell2mat(data(2:end, 2:5)); 

% Average any repeated value
[strains, hit_data] = grpstats(hit_data, strains, {'gname','mean'});

% Prepare final dataset
lis_bobek_2013.orfs = strains;
lis_bobek_2013.ph = strcat(phenotypes, '; ', treatments);
lis_bobek_2013.data = hit_data;

%% Save

save('./lis_bobek_2013.mat','lis_bobek_2013');

fid = fopen('./lis_bobek_2013.txt','w');
write_matrix_file(fid, lis_bobek_2013.orfs, lis_bobek_2013.ph, lis_bobek_2013.data);
fclose(fid);

end
