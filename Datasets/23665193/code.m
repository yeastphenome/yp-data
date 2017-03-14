%% Hirasawa~Shimizu, 2013
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hirasawa_shimizu_2013.pmid = 23665193;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hirasawa_shimizu_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/data_from_S1_PDF.xlsx', 'Mutants');
data.raw(1,:) = [];

% Get the genes from the data
data2.genenames = [data.raw(:,1); data.raw(:,4); data.raw(:,7); data.raw(:,10)];

% Get the data
data2.data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];

% Clean gene names
data2.genenames = clean_genename(data2.genenames);

% Remove all names that are not genes
inds = find(cellfun(@isnumeric, data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

% Remove all names that are not genes
inds = find(~is_genename(data2.genenames) & ~is_orf(data2.genenames));
data2.genenames(inds) = [];
data2.data(inds,:) = [];

% Translate the genes
[data2.orfs, translated] = translate(data2.genenames);
data2.orfs(~translated) = [];
data2.data(~translated, :) = [];

data2.data = cell2mat(data2.data);
data2.data = nanmean(data2.data,2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [94];

%% Load controls

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/data_from_S1_PDF.xlsx', 'CTRL');
data.raw(1,:) = [];
ctrl_data = [data.raw(:,2:3); data.raw(:,5:6); data.raw(:,8:9); data.raw(:,11:12)];
ctrl_data = reshape(ctrl_data,[],1);

inds = find(~cellfun(@isnumeric, ctrl_data));
ctrl_data(inds) = [];
ctrl_data = cell2mat(ctrl_data);
ctrl_data = nanmean(ctrl_data);

data2.data_norm = data2.data ./ ctrl_data;

% Average replicates
[t,t2] = grpstats(data2.data_norm, data2.orfs, {'gname','mean'});

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hirasawa_shimizu_2013.orfs = t;
hirasawa_shimizu_2013.ph = hit_data_names;
hirasawa_shimizu_2013.data = t2;
hirasawa_shimizu_2013.dataset_ids = hit_data_ids;

%% Save

save('./hirasawa_shimizu_2013.mat','hirasawa_shimizu_2013');

%% Print out

fid = fopen('./hirasawa_shimizu_2013.txt','w');
write_matrix_file(fid, hirasawa_shimizu_2013.orfs, hirasawa_shimizu_2013.ph, hirasawa_shimizu_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hirasawa_shimizu_2013)
end
end
