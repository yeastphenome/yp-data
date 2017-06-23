%% Endo~Shima, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
endo_shima_2008.pmid = 18471310;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(endo_shima_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/13068_2007_3_MOESM1_ESM.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(5:end,2);

% Get the data itself
hit_data = data(5:end, 4); 
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Remove empty strain indices
ind = find(cellfun(@isnumeric, hit_strains));
hit_strains(ind) = [];
hit_data(ind) = [];

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% Clean up the data
hit_data = cell2mat(hit_data);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11827];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
endo_shima_2008.orfs = hit_strains;
endo_shima_2008.ph = hit_data_names;
endo_shima_2008.data = hit_data;
endo_shima_2008.dataset_ids = hit_data_ids;

%% Save

save('./endo_shima_2008.mat','endo_shima_2008');

%% Print out

fid = fopen('./endo_shima_2008.txt','w');
write_matrix_file(fid, endo_shima_2008.orfs, endo_shima_2008.ph, endo_shima_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(endo_shima_2008)
end

end