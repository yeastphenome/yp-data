%% Gresham~Botstein, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gresham_botstein_2011.pmid = 20944018;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gresham_botstein_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%
phenotypes = {'Death rate (%/hr)'};
treatments = {'phosphate starvation';'leucine starvation'};

% Dataset #1
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/TABLES4.xlsx', 'phoAbs.txt');

data.genenames = data.raw(2:end, 1);
data.data = cell2mat(data.raw(2:end, 2));
data.genenames_noannot = cell(size(data.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data.genenames)
    C = regexp(data.genenames{i},'_','split');
    data.genenames_noannot{i} = C{1};
end

data.genenames_noannot = clean_genename(data.genenames_noannot);

[data.orfs, translated] = translate(data.genenames_noannot);
data.orfs(~translated) = [];
data.data(~translated,:) = [];


% Dataset #2
[FILENAMES{end+1}, data2.raw] = read_data('xlsread','./raw_data/TABLES5.xlsx', 'leuAbs.txt');

data2.genenames = data2.raw(2:end, 1);
data2.data = cell2mat(data2.raw(2:end, 2));
data2.genenames_noannot = cell(size(data2.genenames));
% Eliminate the "_p" suffix from the genenames
for i = 1 : length(data2.genenames)
    C = regexp(data2.genenames{i},'_','split');
    data2.genenames_noannot{i} = C{1};
end

data2.genenames_noannot = clean_genename(data2.genenames_noannot);

[data2.orfs, translated] = translate(data2.genenames_noannot);
data2.orfs(~translated) = [];
data2.data(~translated,:) = [];


% Average data for identical ORFs that appear multiple times
[data.orfs_u,data.data_avg] = grpstats(data.data, data.orfs, {'gname','mean'});
[data2.orfs_u,data2.data_avg] = grpstats(data2.data, data2.orfs, {'gname','mean'});

hit_orfs = unique([data.orfs_u; data2.orfs_u]);
hit_data = nan(length(hit_orfs),2);

[~,ind1,ind2] = intersect(hit_orfs, data.orfs_u);
hit_data(ind1,1) = data.data_avg(ind2);
[~,ind1,ind2] = intersect(hit_orfs, data2.orfs_u);
hit_data(ind1,2) = data2.data_avg(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [93; 451];

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

gresham_botstein_2011.orfs = hit_orfs;
gresham_botstein_2011.data = hit_data;
gresham_botstein_2011.ph = hit_data_names;
gresham_botstein_2011.dataset_ids = hit_data_ids;

%% Save

save('./gresham_botstein_2011.mat','gresham_botstein_2011');

%% Print out

fid = fopen('./gresham_botstein_2011.txt','w');
write_matrix_file(fid, gresham_botstein_2011.orfs, gresham_botstein_2011.ph, gresham_botstein_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gresham_botstein_2011)
end

end
