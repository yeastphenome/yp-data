%% Gatbonton~Bedalov, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gatbonton_bedalov_2006.pmid = 16552446;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gatbonton_bedalov_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gatbonton_bedalov_2006_hits.xlsx', 'Sheet1');

hit_genenames = clean_genename(data.raw(:,1));

[hits_orfs, translated] = translate(data.raw(:,1));
hits_orfs(~translated) = [];
data.raw(~translated,:) = [];

hits_scores = cell2mat(data.raw(:,3));

% Converting the scores such that:
% a) 1 = weakest phenotype, 3 = strongest phenotype
% b) long telomeres = positive scores, short telomeres = negative scores
hits_scores = abs(hits_scores-4);
inds = find(strcmp('S', data.raw(:,2)));
hits_scores(inds) = -hits_scores(inds);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [104];

%% Load tested genes

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/genelist_altered_020806.xlsx', 'mat alpha copy.txt');

tested_orfs = clean_orf(data.raw(2:end,1));
tested_orfs(strcmp('YYKL138C', tested_orfs)) = {'YKL138C'};

inds = find(~is_orf(tested_orfs));
tested_orfs(inds) = [];

tested_orfs = unique(tested_orfs);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hits_orfs, tested_orfs);
tested_orfs = [tested_orfs; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

gatbonton_bedalov_2006.orfs = tested_orfs;
gatbonton_bedalov_2006.data = zeros(length(tested_orfs), length(hit_data_names));
[~,ind1,ind2] = intersect(gatbonton_bedalov_2006.orfs, hits_orfs);
gatbonton_bedalov_2006.data(ind1,:) = hits_scores(ind2,:);

gatbonton_bedalov_2006.ph = hit_data_names;
gatbonton_bedalov_2006.dataset_ids = hit_data_ids;

%% Save

save('./gatbonton_bedalov_2006.mat','gatbonton_bedalov_2006');

%% Print out

fid = fopen('./gatbonton_bedalov_2006.txt','w');
write_matrix_file(fid, gatbonton_bedalov_2006.orfs, gatbonton_bedalov_2006.ph, gatbonton_bedalov_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gatbonton_bedalov_2006)
end

end
