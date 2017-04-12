%% Burston~Conibear, 2009
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
burston_conibear_2009.pmid = 19506040;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(burston_conibear_2009.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/JCB_200811116_TS1.xlsx', 'TableS1');

hits_orfs = data.raw(6:end,2);
hits_data = data.raw(6:end,4:5);

hits_orfs = clean_orf(hits_orfs);

hits_orfs = translate(hits_orfs);

inds = find(~is_orf(hits_orfs));
disp(hits_orfs(inds))

hits_data = cell2mat(hits_data);

% If the same strain is present more than once, average its values
[hits_orfs, hits_data] = grpstats(hits_data, hits_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [123; 536];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
burston_conibear_2009.orfs = hits_orfs;
burston_conibear_2009.ph = hit_data_names;
burston_conibear_2009.data = hits_data;
burston_conibear_2009.dataset_ids = hit_data_ids;

%% Save

save('./burston_conibear_2009.mat','burston_conibear_2009');

%% Print out

fid = fopen('./burston_conibear_2009.txt','w');
write_matrix_file(fid, burston_conibear_2009.orfs, burston_conibear_2009.ph, burston_conibear_2009.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(burston_conibear_2009)
end

end
