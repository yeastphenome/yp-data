%% Johnson~Wu, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
johnson_wu_2015.pmid = 26688044;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras
% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(johnson_wu_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data
[FILENAMES{end+1}, data] = read_data('xlsread', './raw_data/c5mt00261c1.xlsx', 'Sheet1');

% Get the list of ORFs and the corresponding data
hit_strains = data(4:end, 1);

% Clean up ORFs
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

% Get data from hits
hit_data = cell2mat(data(4:end, 2:7)); 
wt_data = cell2mat(data(3,2:7));

% Standardize the data based on 2 controls: water and wildtype
water = hit_data(:,2)-hit_data(:,1);
wt_water = wt_data(2)-wt_data(1);
water = water ./ wt_water;

concentration_one = hit_data(:,4)-hit_data(:,3);
wt_concentration_one = wt_data(4)-wt_data(3);
concentration_one = concentration_one ./ wt_concentration_one;
concentration_one = concentration_one ./ water;

concentration_two = hit_data(:,6)-hit_data(:,5);
wt_concentration_two = wt_data(6)-wt_data(5);
concentration_two = concentration_two ./ wt_concentration_two;
concentration_two = concentration_two ./ water;

% Add to make final data
hit_data = [concentration_one concentration_two];

% Remove the Inf
hit_data(isinf(hit_data)) = NaN;

% Average any repeated value
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [697; 698];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
johnson_wu_2015.orfs = hit_strains;
johnson_wu_2015.ph = hit_data_names;
johnson_wu_2015.data = hit_data;
johnson_wu_2015.dataset_ids = hit_data_ids;

%% Save

save('./johnson_wu_2015.mat','johnson_wu_2015');

%% Print out

fid = fopen('./johnson_wu_2015.txt','w');
write_matrix_file(fid, johnson_wu_2015.orfs, johnson_wu_2015.ph, johnson_wu_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(johnson_wu_2015)
end

end