%% Garay~DeLuna, 2014
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
garay_deluna_2014.pmid = 24586198;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(garay_deluna_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load hits
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/journal.pgen.1004168.s013.xlsx', 'Genome-wide CLS screen');
hits_orfs = data.raw(2:end,1);
hits_scores = data.raw(2:end,3);

hits_orfs = clean_orf(hits_orfs);
hits_orfs(strcmp('YLR287-A', hits_orfs)) = {'YLR287C-A'};

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds));

hits_scores = cell2mat(hits_scores);

[hits_orfs, hits_scores] = grpstats(hits_scores,hits_orfs,{'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [416];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

garay_deluna_2014.orfs = hits_orfs;
garay_deluna_2014.data = hits_scores;
garay_deluna_2014.ph = hit_data_names;
garay_deluna_2014.dataset_ids = hit_data_ids;

%% Save

save('./garay_deluna_2014.mat','garay_deluna_2014');

%% Print out

fid = fopen('./garay_deluna_2014.txt','w');
write_matrix_file(fid, garay_deluna_2014.orfs, garay_deluna_2014.ph, garay_deluna_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(garay_deluna_2014)
end

end
