%% Zhao~Jiang, 2013

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhao_jiang_2013.pmid = 23026396;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zhao_jiang_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hit Data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
hit_strains = data(3:end,1);

% Get the data itself
hit_data = -ones(length(hit_strains), 1);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1309];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/DELETION LIBRARY.xlsx', 'DELETION LIBRARY');

% Get the list of ORFs and the correponding data 
tested_strains = data(3:end,2);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
tested_strains(ismember(tested_strains, {'YELOO1C'})) = {'YEL001C'};
inds = find(~is_orf(tested_strains));
tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
zhao_jiang_2013.orfs = tested_strains;
zhao_jiang_2013.ph = hit_data_names;
zhao_jiang_2013.data = zeros(length(zhao_jiang_2013.orfs),length(zhao_jiang_2013.ph));
zhao_jiang_2013.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, zhao_jiang_2013.orfs);
zhao_jiang_2013.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./zhao_jiang_2013.mat','zhao_jiang_2013');

%% Print out

fid = fopen('./zhao_jiang_2013.txt','w');
write_matrix_file(fid, zhao_jiang_2013.orfs, zhao_jiang_2013.ph, zhao_jiang_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zhao_jiang_2013)
end

end

