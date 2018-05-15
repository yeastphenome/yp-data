%% Hartman~Tippery, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hartman_tippery_2004.pmid = 15239834;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hartman_tippery_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%%

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/gb-2004-5-7-r49-s7.xlsx', 'data');

ind_orf = strmatch('ORF', data.raw(1,:));
hit_orfs = data.raw(2:end, ind_orf);

ind_data1 = strmatch('No HU- Growth Index', data.raw(1,:));
ind_data2 = strmatch('50mM HU Growth Index', data.raw(1,:));
ind_data3 = strmatch('150 mM HU Growth Index', data.raw(1,:));

hit_data = data.raw(2:end, [ind_data1 ind_data2 ind_data3]);

% Unfortunate case: ~70 ORFs have an extra symbol ("b", "c", etc.) that
% looks like it is part of the ORF, but it isn't. 
inds = find(cellfun(@length, hit_orfs) == 8);   % ~70 ORFs have a "-" missing
for i = 1 : length(inds)
    hit_orfs{inds(i)} = hit_orfs{inds(i)}(1:7);
end

% Eliminate anything that doesn't look like an ORF
hit_orfs = clean_orf(hit_orfs);

hit_orfs(strcmp('YOR298C-AB', hit_orfs)) = {'YOR298C-A'};

inds = find(cellfun(@isnumeric, hit_orfs));
hit_orfs(inds) = [];
hit_data(inds,:) = [];

hit_orfs = translate(hit_orfs);

inds = find(~is_orf(hit_orfs));
disp(hit_orfs(inds));

hit_orfs(inds) = [];
hit_data(inds,:) = [];

% Make sure all the data are numbers
inds = find(~cellfun(@isnumeric, hit_data));
hit_data(inds) = {NaN};
hit_data = cell2mat(hit_data);

% Normalize by UNTREATED sample
hit_data(:,2) = hit_data(:,2) - hit_data(:,1);
hit_data(:,3) = hit_data(:,3) - hit_data(:,1);

% Average data for identical ORFs that appear multiple times
[hit_orfs,hit_data] = grpstats(hit_data, hit_orfs, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16186; 52; 53];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

hartman_tippery_2004.orfs = hit_orfs;
hartman_tippery_2004.data = hit_data;
hartman_tippery_2004.ph = hit_data_names;
hartman_tippery_2004.dataset_ids = hit_data_ids;

%% Save

save('./hartman_tippery_2004.mat','hartman_tippery_2004');

%% Print out

fid = fopen('./hartman_tippery_2004.txt','w');
write_matrix_file(fid, hartman_tippery_2004.orfs, hartman_tippery_2004.ph, hartman_tippery_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hartman_tippery_2004)
end


end
