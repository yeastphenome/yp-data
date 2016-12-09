%% Banuelos~Gharakhanian, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
banuelos_gharakhanian_2010.pmid = 20043226;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(banuelos_gharakhanian_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,2);
hit_data = cell2mat(hit_data);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

hit_strains(strcmp('YM1010W-A', hit_strains)) = {'YML010W-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [154];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
banuelos_gharakhanian_2010.orfs = hit_strains;
banuelos_gharakhanian_2010.ph = hit_data_names;
banuelos_gharakhanian_2010.data = hit_data;
banuelos_gharakhanian_2010.dataset_ids = hit_data_ids;

%% Save
save('./banuelos_gharakhanian_2010.mat','banuelos_gharakhanian_2010');

%% Print out

fid = fopen('./banuelos_gharakhanian_2010.txt','w');
write_matrix_file(fid, banuelos_gharakhanian_2010.orfs, banuelos_gharakhanian_2010.ph, banuelos_gharakhanian_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(banuelos_gharakhanian_2010)
end

end

