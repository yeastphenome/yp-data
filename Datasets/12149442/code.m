%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
hanway_romesberg_2002.pmid = 12149442;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(hanway_romesberg_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data (UV)

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/uv_hits.xlsx', 'UV');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);

% Get the data itself
hit_data = data(:,2:9);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'NTG1D'})) = {'NTG1'};

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cell2mat(hit_data)-1;   % Coverting to scale where 0 = wt, and negative values = sensitivity

%% Load data (MMS)

[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/uv_hits.xlsx', 'MMS');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2(:,1);

% Get the data itself
hit_data2 = data2(:,2:3);
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_genename(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains2(inds) = [];
hit_data2(inds,:) = [];

hit_data2 = cell2mat(hit_data2)-1;   % Coverting to scale where 0 = wt, and negative values = sensitivity

%%

hit_strains_all = unique([hit_strains; hit_strains2]);
hit_data_all = nan(length(hit_strains_all), 10);

[~,ind1,ind2] = intersect(hit_strains, hit_strains_all);
hit_data_all(ind2,1:8) = hit_data(ind1,:);

[~,ind1,ind2] = intersect(hit_strains2, hit_strains_all);
hit_data_all(ind2,9:10) = hit_data2(ind1,:);

%% 
% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data_all, hit_strains_all, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [474 11850 11851 11852 11853 11854 11855 11856 4951 4967]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
hanway_romesberg_2002.orfs = hit_strains;
hanway_romesberg_2002.ph = hit_data_names;
hanway_romesberg_2002.data = hit_data;
hanway_romesberg_2002.dataset_ids = hit_data_ids;

%% Save

save('./hanway_romesberg_2002.mat','hanway_romesberg_2002');

%% Print out

fid = fopen('./hanway_romesberg_2002.txt','w');
write_matrix_file(fid, hanway_romesberg_2002.orfs, hanway_romesberg_2002.ph, hanway_romesberg_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(hanway_romesberg_2002)
end

end

