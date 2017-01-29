%% Wang~Zhou, 2007
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wang_zhou_2007.pmid = 17308930;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(wang_zhou_2007.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data: EDTA data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mutants sensitive to metal scarcity.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
edta_strains = data(2:end,1);

% Get the data itself
indx = find(~cellfun(@isempty, regexp(data(2:end,3), 'HS', 'match')));
edta_data = zeros(length(edta_strains), 1)-1;
edta_data(indx) = -2;

% Eliminate all white spaces & capitalize
edta_strains = clean_orf(edta_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(edta_strains));
edta_strains(inds) = [];
edta_data(inds) = [];

%% Load the data: metal data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mutants sensitive to metal excess.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
metal_strains = data(2:end,1);

% Get the data itself
metals = data(2:end, 3:6);
metal_data = zeros(size(metals));
for i = 1:size(metals, 2)
    sample = metals(:,i);
    C = cellfun(@(x) x(1), sample, 'un', 0);
    indx1 = find(~cellfun(@isempty, regexp(sample, 'HS', 'match')));
    metal_data(indx1,i) = -2;
    indx2 = find(~cellfun(@isempty, regexp(C, 'S', 'match')));
    metal_data(indx2,i) = -1;
end

% Eliminate all white spaces & capitalize
metal_strains = clean_orf(metal_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(metal_strains));
metal_strains(inds) = [];
metal_data(inds, :) = [];

%% Combine the data
hit_strains = unique([metal_strains; edta_strains]);
hit_data = nan(length(hit_strains), 5);
[~, ind1, ind2] = intersect(hit_strains,edta_strains);
hit_data(ind1, 1) = edta_data(ind2);
[~, ind1, ind2] = intersect(hit_strains,metal_strains);
hit_data(ind1, 2:5) = metal_data(ind2, :);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5246; 5247; 5248; 5249; 5250];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
wang_zhou_2007.orfs = hit_strains;
wang_zhou_2007.ph = hit_data_names;
wang_zhou_2007.data = hit_data;
wang_zhou_2007.dataset_ids = hit_data_ids;

%% Save

save('./wang_zhou_2007.mat','wang_zhou_2007');

%% Print out

fid = fopen('./wang_zhou_2007.txt','w');
write_matrix_file(fid, wang_zhou_2007.orfs, wang_zhou_2007.ph, wang_zhou_2007.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(wang_zhou_2007)
end

end
