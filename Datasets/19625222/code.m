%% de graaf~McCullough, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
de_graaf_mccullough_2009.pmid = 19625222;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(de_graaf_mccullough_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,2);

% Split the strings
C = cellfun(@(s) strsplit(s, ' '), hit_strains, 'UniformOutput', 0);
hit_strains = [C{1}'; C{2}'; C{3}'; C{4}'];

% Get the data itself
hit_data = zeros(size(hit_strains))-1;

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [161];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
de_graaf_mccullough_2009.orfs = hit_strains;
de_graaf_mccullough_2009.ph = hit_data_names;
de_graaf_mccullough_2009.data = hit_data;
de_graaf_mccullough_2009.dataset_ids = hit_data_ids;

%% Save

save('./de_graaf_mccullough_2009.mat','de_graaf_mccullough_2009');

%% Print out

fid = fopen('./de_graaf_mccullough_2009.txt','w');
write_matrix_file(fid, de_graaf_mccullough_2009.orfs, de_graaf_mccullough_2009.ph, de_graaf_mccullough_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(de_graaf_mccullough_2009)
end

end

