%% Czyz~Zaremberg, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
czyz_zaremberg_2013.pmid = 23344949;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(czyz_zaremberg_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
hit_strains = data(2:end,2);

% Get the data itself
hit_data = data(2:end,1);
   
% Eliminate all white spaces & capitalize
hit_strains = regexprep(hit_strains,'\(\w+\)', '');
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% Fix up the data
hit_data = strrep(hit_data, '+++', '-3');
hit_data = strrep(hit_data, '++', '-2');
hit_data = strrep(hit_data, '+', '-1');
hit_data = str2double(hit_data);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [774];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
czyz_zaremberg_2013.orfs = hit_strains;
czyz_zaremberg_2013.ph = hit_data_names;
czyz_zaremberg_2013.data = hit_data;
czyz_zaremberg_2013.dataset_ids = hit_data_ids;

%% Save

save('./czyz_zaremberg_2013.mat','czyz_zaremberg_2013');

%% Print out

fid = fopen('./czyz_zaremberg_2013.txt','w');
write_matrix_file(fid, czyz_zaremberg_2013.orfs, czyz_zaremberg_2013.ph, czyz_zaremberg_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(czyz_zaremberg_2013)
end

end