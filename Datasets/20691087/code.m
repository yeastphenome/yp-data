%% Alamgir~Golshani, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alamgir_golshani_2010.pmid = 20691087;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(alamgir_golshani_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%
[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/1472-6769-10-6-s1.xlsx', 'Raw genome-wide data');

% Eliminate anything that doesn't look like an ORF
hit_orfs = data.raw(:,1);
hit_data = data.raw(:,3:14);

hit_orfs = clean_orf(hit_orfs);

hit_orfs(strcmp('YPL072WA', hit_orfs)) = {'YPL072W-A'};
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

hit_orfs(inds) = [];
hit_data(inds,:) = []; 

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% From looking at the data, it appears that the numbers indicate
% sensitivity, not growth, so we will revert their sign, such that higher
% values indicate better growth and lower values indicate the opposite.
hit_data = -hit_data;

% Average data for identical ORFs that appear multiple times
[hit_orfs,hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [9 9 9 10 10 10 7 7 7 8 8 8];
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare the final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

alamgir_golshani_2010.orfs = hit_orfs;
alamgir_golshani_2010.data = hit_data;
alamgir_golshani_2010.ph = hit_data_names;
alamgir_golshani_2010.dataset_ids = hit_data_ids;

%% Save

save('./alamgir_golshani_2010.mat','alamgir_golshani_2010');

%% Print out

fid = fopen('./alamgir_golshani_2010.txt','w');
write_matrix_file(fid, alamgir_golshani_2010.orfs, alamgir_golshani_2010.ph, alamgir_golshani_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(alamgir_golshani_2010)
end

end
