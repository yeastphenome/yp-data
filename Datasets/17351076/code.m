%% Ohnuki~Ohya, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ohnuki_ohya_2007.pmid = 17351076;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ohnuki_ohya_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table2.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,2);

% Get the data itself
hit_data = data(2:end, 3); 

% Split the strains and get the ORF
C = regexp(hit_strains, '/', 'split');
hit_strains = {};
for i = 1:length(C)
    if length(C{i}) == 2
        hit_strains(end+1,:) = C{i}(2);
    else
        hit_strains(end+1,:) = C{i}(1);
    end
end
    
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YOR331CF'})) = {'YOR331C'};
hit_strains(ismember(hit_strains, {'YPR099CG'})) = {'YPR099C'};

% "D" and "E" are subscript annotations here, not ORF name variants
hit_strains(ismember(hit_strains, {'YEL045C-D'})) = {'YEL045C'};
hit_strains(ismember(hit_strains, {'YKL118W-E'})) = {'YKL118W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];  
hit_data(inds, :) = [];

% Transform the data
final_hit_data = zeros(length(hit_data),1);
letters = {'C','B','A'};
for i = 1:length(letters)
    ind = strcmp(letters{i}, hit_data);
    final_hit_data(ind) = -i;
end

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11826];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ohnuki_ohya_2007.orfs = hit_strains;
ohnuki_ohya_2007.ph = hit_data_names;
ohnuki_ohya_2007.data = final_hit_data;
ohnuki_ohya_2007.dataset_ids = hit_data_ids;

%% Save

save('./ohnuki_ohya_2007.mat','ohnuki_ohya_2007');

%% Print out

fid = fopen('./ohnuki_ohya_2007.txt','w');
write_matrix_file(fid, ohnuki_ohya_2007.orfs, ohnuki_ohya_2007.ph, ohnuki_ohya_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ohnuki_ohya_2007)
end

end