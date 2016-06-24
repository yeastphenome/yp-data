%% Chan~Zheng, 2000
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
chan_zheng_2000.pmid = 11078525;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(chan_zheng_2000.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/chan_zheng_2000_HAP.xlsx', 'Sheet1');

hit_orfs = data.raw(:,1);
hit_data = cell2mat(data.raw(:,2));
hit_data = log10(hit_data);

hit_orfs = clean_orf(hit_orfs);

inds = find(~is_orf(hit_orfs));
disp(inds);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [1];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

chan_zheng_2000.orfs = hit_orfs;
chan_zheng_2000.data = hit_data;
chan_zheng_2000.ph = hit_data_names;
chan_zheng_2000.dataset_ids = hit_data_ids;

%% Save

save('./chan_zheng_2000.mat','chan_zheng_2000');

%% Print out

fid = fopen('./chan_zheng_2000.txt','w');
write_matrix_file(fid, chan_zheng_2000.orfs, chan_zheng_2000.ph, chan_zheng_2000.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(chan_zheng_2000)
end