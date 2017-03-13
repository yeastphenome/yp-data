%% Stepchenkova~Pavlov, 2005
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
stepchenkova_pavlov_2005.pmid = 15932646;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(stepchenkova_pavlov_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/stepchenkova_pavlov_2005_hits.xlsx', 'Sheet1');

hits_genenames = data.raw(:,1);
hits_scores = data.raw(:,2:3);

hits_genenames = clean_genename(hits_genenames);

inds = find(cellfun(@isnumeric, hits_genenames));
hits_genenames(inds) = [];
hits_scores(inds,:) = [];

hits_orfs = translate(hits_genenames);

hits_scores = cell2mat(hits_scores);

% Normalize to WT
inds = find(strcmpi('WT', hits_orfs));
hits_scores = hits_scores ./ repmat(hits_scores(inds,:), length(hits_orfs),1);
hits_scores(inds,:) = [];
hits_orfs(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [186 423]';

%% Load tested genes

[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/mat_alpha_041902.txt', '%*s %s %*s %*s %*s %*s %*s %*s', 'delimiter', '\t');

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds)); 

tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

tested_orfs = [tested_orfs; missing];   % 1 missing strain added

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
stepchenkova_pavlov_2005.orfs = tested_orfs;
stepchenkova_pavlov_2005.ph = hit_data_names;
stepchenkova_pavlov_2005.data = zeros(length(stepchenkova_pavlov_2005.orfs),length(stepchenkova_pavlov_2005.ph));
stepchenkova_pavlov_2005.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, stepchenkova_pavlov_2005.orfs);
stepchenkova_pavlov_2005.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./stepchenkova_pavlov_2005.mat','stepchenkova_pavlov_2005');

%% Print out

fid = fopen('./stepchenkova_pavlov_2005.txt','w');
write_matrix_file(fid, stepchenkova_pavlov_2005.orfs, stepchenkova_pavlov_2005.ph, stepchenkova_pavlov_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(stepchenkova_pavlov_2005)
end

end
