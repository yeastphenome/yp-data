%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
sousa_sousa_2013.pmid = 24286259;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(sousa_sousa_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/12864_2013_5541_MOESM3_ESM.xlsx', 'Folha2');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,2);
hit_data = data(:,3); 
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds,:) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YDL230'})) = {'YDL230W'};
hit_strains(ismember(hit_strains, {'YDL178'})) = {'YDL178W'};
hit_strains(ismember(hit_strains, {'YLR124W-'})) = {'YLR124W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data(find(strcmp(hit_data, '+'))) = {1};
hit_data(find(strcmp(hit_data, '-'))) = {-1};
hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data, count] = grpstats(hit_data, hit_strains, {'gname','mean', 'numel'});


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11812];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
sousa_sousa_2013.orfs = hit_strains;
sousa_sousa_2013.ph = hit_data_names;
sousa_sousa_2013.data = hit_data;
sousa_sousa_2013.dataset_ids = hit_data_ids;

%% Save

save('./sousa_sousa_2013.mat','sousa_sousa_2013');

%% Print out

fid = fopen('./sousa_sousa_2013.txt','w');
write_matrix_file(fid, sousa_sousa_2013.orfs, sousa_sousa_2013.ph, sousa_sousa_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(sousa_sousa_2013)
end

end

