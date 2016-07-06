%% Blackburn~Avery, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available

blackburn_avery_2003.pmid = 12543677;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(blackburn_avery_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/blackburn_avery_2003_data.xlsx', 'data.txt');

treatments = data.raw(1,2:8)';
data.raw(1,:) = [];

% Eliminate white spaces before/after ORF
hit_orfs = data.raw(:,1);
hit_data = data.raw(:,2:8);

hit_orfs = clean_orf(hit_orfs);

% Eliminate everything that doesn't look like an ORF
inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

% Replace 'Inf' with Inf
hit_data(~cellfun(@isnumeric, hit_data)) = {5120};  % MIC=Inf set to 10X the maximum MIC detected
hit_data = cell2mat(hit_data);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [391; 393; 395; 69; 392; 394; 396];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

blackburn_avery_2003.orfs = hit_orfs;
blackburn_avery_2003.data = hit_data;
blackburn_avery_2003.ph = hit_data_names;
blackburn_avery_2003.dataset_ids = hit_data_ids;

%% Save

save('./blackburn_avery_2003.mat','blackburn_avery_2003');

%% Print out

fid = fopen('./blackburn_avery_2003.txt','w');
write_matrix_file(fid, blackburn_avery_2003.orfs, blackburn_avery_2003.ph, blackburn_avery_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(blackburn_avery_2003)
end

end

