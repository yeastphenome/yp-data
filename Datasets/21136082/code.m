%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
villa_garcia_henry_2011.pmid = 21136082;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(villa_garcia_henry_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/NIHMS265254-supplement-2.xlsx', 'Table 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,3:5);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);
inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data = cellfun(@upper, hit_data, 'UniformOutput', 0);

% Convert the phenotypes
hit_data(find(strcmp('+', hit_data))) = {0};
hit_data(find(strcmp('VW', hit_data))) = {-1};
hit_data(find(strcmp('W', hit_data))) = {-2};
hit_data(find(strcmp('S', hit_data))) = {-3};
hit_data(find(strcmp('NS', hit_data))) = {NaN};
hit_data(find(strcmp('NG', hit_data))) = {NaN};

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16191; 16192; 16193];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
villa_garcia_henry_2011.orfs = hit_strains;
villa_garcia_henry_2011.ph = hit_data_names;
villa_garcia_henry_2011.data = hit_data;
villa_garcia_henry_2011.dataset_ids = hit_data_ids;

%% Save

save('./villa_garcia_henry_2011.mat','villa_garcia_henry_2011');

%% Print out

fid = fopen('./villa_garcia_henry_2011.txt','w');
write_matrix_file(fid, villa_garcia_henry_2011.orfs, villa_garcia_henry_2011.ph, villa_garcia_henry_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(villa_garcia_henry_2011)
end

end

