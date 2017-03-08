%% Reiner~Schneiter, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
reiner_schneiter_2006.pmid = 16251356;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(reiner_schneiter_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested genes

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/BY4741-MATa COLLECTION.xls', 'chr11_1yes');

tested_orfs = data.raw(2:end,2);

tested_orfs = clean_orf(tested_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_orfs));
disp(tested_orfs(inds));  

tested_orfs(inds) = [];
tested_orfs = unique(tested_orfs);

%% Load hits

[FILENAMES{end+1}, hits] = read_data('textread','./raw_data/reiner_scheiter_2006_hits.txt', '%s');

% This list of ORFs is lacking the last character (published that way), so
% we have to match it to the list of tested strains.
for i = 1 : length(hits)
    inds = find(strncmp(hits{i}, tested_orfs, length(hits{i})));
    if length(inds) == 1
        hits_orfs(i) = tested_orfs(inds);
    else
        fprintf('%s\t%d\n', hits{i}, length(inds));
    end
end

% Two ORFs (YBR039W and YNL243W) could not be found in the list of tested
% strains. So we have to remove them.
hits_orfs(ismember(hits, {'YBR039','YNL243'})) = [];

hits_orfs = unique(hits_orfs);

% Assign a -1 to all the strains.
hits_scores = zeros(length(hits_orfs),1)-1;

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [175];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
reiner_schneiter_2006.orfs = tested_orfs;
reiner_schneiter_2006.ph = hit_data_names;
reiner_schneiter_2006.data = zeros(length(reiner_schneiter_2006.orfs),length(reiner_schneiter_2006.ph));
reiner_schneiter_2006.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, reiner_schneiter_2006.orfs);
reiner_schneiter_2006.data(ind2,:) = hits_scores(ind1,:);

%% Save

save('./reiner_schneiter_2006.mat','reiner_schneiter_2006');

%% Print out

fid = fopen('./reiner_schneiter_2006.txt','w');
write_matrix_file(fid, reiner_schneiter_2006.orfs, reiner_schneiter_2006.ph, reiner_schneiter_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(reiner_schneiter_2006)
end

end
