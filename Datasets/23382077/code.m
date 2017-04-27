%% Huang~Paulovich, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
huang_paulovich_2013.pmid = 23382077;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(huang_paulovich_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

% Load plate maps
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/Mat_a_obs_v4 0.xlsx', 'DATA');
tested_orfs = tested.raw(2:end,2);

% Eliminate all white spaces & capitalize
tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(cellfun(@isnumeric,tested_orfs));
tested_orfs(inds)= [];  
tested_orfs(ismember(tested_orfs, {'YLR287-A'})) = {'YLR287C-A'};

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [145];

%% Load tested data
[FILENAMES{end+1}, hits_genenames] = read_data('textread','./raw_data/huang_paulovich_2013_hits.txt', '%s');

hits_genenames = clean_genename(hits_genenames);
hits_orfs = translate(hits_genenames);

hits_scores = -ones(length(hits_orfs),1);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
huang_paulovich_2013.orfs = tested_orfs;
huang_paulovich_2013.ph = hit_data_names;
huang_paulovich_2013.data = zeros(length(huang_paulovich_2013.orfs),length(huang_paulovich_2013.ph));
huang_paulovich_2013.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, huang_paulovich_2013.orfs);
huang_paulovich_2013.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./huang_paulovich_2013.mat','huang_paulovich_2013');

%% Print out

fid = fopen('./huang_paulovich_2013.txt','w');
write_matrix_file(fid, huang_paulovich_2013.orfs, huang_paulovich_2013.ph, huang_paulovich_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(huang_paulovich_2013)
end

end
