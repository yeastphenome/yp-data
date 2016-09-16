%% McCormick~Kennedy, 2015

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mccormick_kennedy_2015.pmid = 26456335;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mccormick_kennedy_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hit Data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc3.xlsx', 'Table S2');

% Get the list of ORFs and the correponding data 
hit_strains = data(4:end,1);

% Get the data itself
hit_data = cell2mat(data(4:end, 12));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
hit_strains(ismember(hit_strains, {'FMP42'})) = {'YMR221C'};
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [696];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
tested_strains = data(3:end,1);

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% Find anything that doesn't look like an ORF
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
mccormick_kennedy_2015.orfs = tested_strains;
mccormick_kennedy_2015.ph = hit_data_names;
mccormick_kennedy_2015.data = zeros(length(mccormick_kennedy_2015.orfs),length(mccormick_kennedy_2015.ph));
mccormick_kennedy_2015.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, mccormick_kennedy_2015.orfs);
mccormick_kennedy_2015.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./mccormick_kennedy_2015.mat','mccormick_kennedy_2015');

%% Print out

fid = fopen('./mccormick_kennedy_2015.txt','w');
write_matrix_file(fid, mccormick_kennedy_2015.orfs, mccormick_kennedy_2015.ph, mccormick_kennedy_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mccormick_kennedy_2015)
end

end

