%% Zhou~Costa, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zhou_costa_2009.pmid = 19631266;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zhou_costa_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load hits

[FILENAMES{end+1}, hits.raw] = read_data('xlsread','./raw_data/table.xlsx', 'table.csv');
hits_orfs = hits.raw(4:end,1);
hits_scores = hits.raw(4:end,4);

hits_orfs = clean_orf(hits_orfs);

hits_orfs = translate(hits_orfs);

hits_orfs(strcmp('YAR002AW', hits_orfs)) = {'YAR002W-A'};

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds))

hits_scores(strcmp('Above 5', hits_scores)) = {5};
hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

wt = 4.47;  %IC50 from paper
hits_scores = hits_scores - wt; % negative scores = more sensitive than WT, and viceversa.

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [115];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zhou_costa_2009.orfs = hits_orfs;
zhou_costa_2009.ph = hit_data_names;
zhou_costa_2009.data = hits_scores;
zhou_costa_2009.dataset_ids = hit_data_ids;

%% Save

save('./zhou_costa_2009.mat','zhou_costa_2009');

%% Print out

fid = fopen('./zhou_costa_2009.txt','w');
write_matrix_file(fid, zhou_costa_2009.orfs, zhou_costa_2009.ph, zhou_costa_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zhou_costa_2009)
end

end
