%% Yuen~Spencer, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yuen_spencer_2007.pmid = 17360454;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yuen_spencer_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/10642Tables_3-9.xlsx', 'SI Table 3');  % hits identified in more than 1 screen
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/10642Tables_3-9.xlsx', 'SI Table 4');  % hits identified in 1 screen only

data = [data1(:,1:6); data2(:,1:6)];

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,2);

% Get the data itself
hit_data = data(2:end,4:6);
hit_data(~cellfun(@isnumeric, hit_data)) = {NaN};
hit_data = cell2mat(hit_data);
hit_data = hit_data - 1;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
[hit_strains, translated, ambiguous] = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4931; 4932; 4933];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yuen_spencer_2007.orfs = hit_strains;
yuen_spencer_2007.ph = hit_data_names;
yuen_spencer_2007.data = hit_data;
yuen_spencer_2007.dataset_ids = hit_data_ids;

%% Save

save('./yuen_spencer_2007.mat','yuen_spencer_2007');

%% Print out

fid = fopen('./yuen_spencer_2007.txt','w');
write_matrix_file(fid, yuen_spencer_2007.orfs, yuen_spencer_2007.ph, yuen_spencer_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yuen_spencer_2007)
end

end