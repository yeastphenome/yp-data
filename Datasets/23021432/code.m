%% Orij~Smits, 2012
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
orij_smits_2012.pmid = 23021432;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(orij_smits_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
% Load hits
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/pH Screen raw.xlsx', 'initial screens');

% Get the list of ORFs and corresponding data
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3:4);

% Eliminate all white spaces & capitalize
hits_orfs = clean_orf(hits_orfs);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hits_orfs));
hits_orfs(inds) = [];
hits_scores(inds, :) = [];

% Translate to correct for old aliases that look like ORFs
[hits_orfs, translated, ambiguous] = translate(hits_orfs);
disp(hits_orfs(find(translated)));
disp(hits_orfs(find(ambiguous)));

% Convert to numeric from cell
hits_scores = cell2mat(hits_scores);
hits_scores = nanmean(hits_scores, 2);

% Get averages
[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [126];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
orij_smits_2012.orfs = hits_orfs;
orij_smits_2012.ph = hit_data_names;
orij_smits_2012.data = hits_scores;
orij_smits_2012.dataset_ids = hit_data_ids;

%% Save

save('./orij_smits_2012.mat','orij_smits_2012');

%% Print out

fid = fopen('./orij_smits_2012.txt','w');
write_matrix_file(fid, orij_smits_2012.orfs, orij_smits_2012.ph, orij_smits_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(orij_smits_2012)
end

end
