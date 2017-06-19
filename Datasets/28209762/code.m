%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
bae_swanson_2017.pmid = 28209762;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(bae_swanson_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);
  
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];

hit_strains = unique(hit_strains);
hit_data = -ones(size(hit_strains));

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11809];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('xlsread','./raw_data/Heterozygous_diploid_obs_v7.0.xlsx', 'Heterozygous Diploid_obs');

tested_strains = tested_strains(2:end,2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YCLO51W'})) = {'YCL051W'};
tested_strains(ismember(tested_strains, {'YHR139C-'})) = {'YHR139C-A'};
tested_strains(ismember(tested_strains, {'YGR122C-'})) = {'YGR122C-A'};

inds = find(~cellfun(@isempty, regexp(tested_strains, '.*1$')));
for i = 1 : length(tested_strains(inds))
    tested_strains{inds(i)} = tested_strains{inds(i)}(1:length(tested_strains{inds(i)})-1);
end

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
bae_swanson_2017.orfs = tested_strains;
bae_swanson_2017.ph = hit_data_names;
bae_swanson_2017.data = zeros(length(bae_swanson_2017.orfs),length(bae_swanson_2017.ph));
bae_swanson_2017.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, bae_swanson_2017.orfs);
bae_swanson_2017.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./bae_swanson_2017.mat','bae_swanson_2017');

%% Print out

fid = fopen('./bae_swanson_2017.txt','w');
write_matrix_file(fid, bae_swanson_2017.orfs, bae_swanson_2017.ph, bae_swanson_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(bae_swanson_2017)
end

end

