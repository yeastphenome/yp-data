%% Yoshikawa~Shimizu, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yoshikawa_shimizu_2009.pmid = 19054128;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yoshikawa_shimizu_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

phenotypes = {'Growth, exponential growth rate'};
treatments = {'UNT';'EtOH, 5%';'EtOH, 8%';'NaCl, 1 M'};

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/FYR_456_sm_tableS1.xlsx', 'data');

hit_orfs = data.raw(3:end,1);

ind_data = [3 4 6 7 9 10 12 13];
hit_data = data.raw(3:end,ind_data);

% Eliminate anything that doesn't look like an ORF
hit_orfs = clean_orf(hit_orfs);
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Average data for identical ORFs that appear multiple times
[hit_orfs,hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [100 100 523 523 4 4 5 5];   % setting UNT to 100 for now
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

% Normalize by UNTREATED sample
ind = find(hit_data_ids == 100);
hit_data = hit_data ./ repmat(hit_data(:,ind),1,4);
hit_data(:,ind) = [];
hit_data_ids(ind) = [];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% Average data for identical ORFs that appear multiple times
yoshikawa_shimizu_2009.orfs = hit_orfs;
yoshikawa_shimizu_2009.data = hit_data;
yoshikawa_shimizu_2009.ph = hit_data_names;
yoshikawa_shimizu_2009.dataset_ids = hit_data_ids;

%% Save

save('./yoshikawa_shimizu_2009.mat','yoshikawa_shimizu_2009');

%% Print out

fid = fopen('./yoshikawa_shimizu_2009.txt','w');
write_matrix_file(fid, yoshikawa_shimizu_2009.orfs, yoshikawa_shimizu_2009.ph, yoshikawa_shimizu_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yoshikawa_shimizu_2009)
end

end
