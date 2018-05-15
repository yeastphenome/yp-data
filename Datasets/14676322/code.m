%% Warringer~Blomberg, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

warringer_blomberg_2003.pmid = 14676322;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(warringer_blomberg_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/LPI NaCl.xlsx', 'LPI');

hit_strains1 = data.raw(5:end,1);
hit_data1 = cell2mat(data.raw(5:end,2:4));

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

hit_strains1(inds) = [];
hit_data1(inds,:) = [];

% Average data for identical ORFs that appear multiple times
[hit_strains1,hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids1 = [50 49 51]';

%% Load control data

[FILENAMES{end+1}, data2.raw] = read_data('xlsread','./raw_data/LSC Reference.xlsx', 'LSC');

hit_strains2 = data2.raw(5:end,1);
hit_data2 = data2.raw(5:end,[4,8,12]);

% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains2(inds) = [];
hit_data2(inds,:) = [];

inds = find(~cellfun(@isnumeric, hit_data2));
hit_data2(inds) = {NaN};

hit_data2 = cell2mat(hit_data2);

% Average data for identical ORFs that appear multiple times
[hit_strains2,hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [16184 16183 16185]';


%% Prepare final dataset

hit_data_ids = [hit_data_ids1; hit_data_ids2];
hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains), length(hit_data_ids));

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1:3) = hit_data1(ind2,:);
[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,4:6) = hit_data2(ind2,:);

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
warringer_blomberg_2003.orfs = hit_strains;
warringer_blomberg_2003.ph = hit_data_names;
warringer_blomberg_2003.data = hit_data;
warringer_blomberg_2003.dataset_ids = hit_data_ids;

%% Save

save('./warringer_blomberg_2003.mat','warringer_blomberg_2003');

%% Print out

fid = fopen('./warringer_blomberg_2003.txt','w');
write_matrix_file(fid, warringer_blomberg_2003.orfs, warringer_blomberg_2003.ph, warringer_blomberg_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(warringer_blomberg_2003)
end

end

