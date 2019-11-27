%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
puddu_jackson_2019.pmid = 31511699;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(puddu_jackson_2019.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/41586_2019_1549_MOESM3_ESM.xlsx', 'SupplementaryTable3');
data = data(35:end,:);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains_ids = data(:,13);
hit_data = data(:,[2:11,14:29,31]);

% Extract gene names
str = 'Del(\d)*_([\w,])*';

hit_strains = cell(length(hit_strains_ids),1);
for i = 1:length(hit_strains_ids)
    [tokens, matches] = regexp(hit_strains_ids{i}, str, 'tokens','match');
    hit_strains{i} = tokens{1}{2};
end

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains); 

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'FLO8'})) = {'YER109C'};
hit_strains(ismember(hit_strains, {'AAD6'})) = {'YFL056C'};
hit_strains(ismember(hit_strains, {'SDL1'})) = {'YIL167W'};
hit_strains(ismember(hit_strains, {'HXT12'})) = {'YIL170W'};
hit_strains(ismember(hit_strains, {'SDC25'})) = {'YLL016W'};
hit_strains(ismember(hit_strains, {'CRS5'})) = {'YOR031W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

% Remove the WTs (for now)
hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in orhit_der)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16322:16348]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
puddu_jackson_2019.orfs = hit_strains;
puddu_jackson_2019.ph = hit_data_names;
puddu_jackson_2019.data = hit_data;
puddu_jackson_2019.dataset_ids = hit_data_ids;

%% Save

save('./puddu_jackson_2019.mat','puddu_jackson_2019');

%% Print out

fid = fopen('./puddu_jackson_2019.txt','w');
write_matrix_file(fid, puddu_jackson_2019.orfs, puddu_jackson_2019.ph, puddu_jackson_2019.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(puddu_jackson_2019)
end

end

