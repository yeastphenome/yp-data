%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
zewaii_huang_2003.pmid = 12615994;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(zewaii_huang_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/0118Table3.xlsx', 'Sheet1');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/0118Table4.xlsx', 'Sheet1');
[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/0118Table5.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = [data1(:,3); data2(:,3); data3(:,3)];

% Get the data itself
hit_data = [data1(:,2); data2(:,2); data3(:,2)];

hit_collection = [data1(:,6); data2(:,6); data3(:,6)];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(strcmp(hit_strains, 'YLR287-A')) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains(inds) = [];
hit_data(inds) = [];
hit_collection(inds) = [];

inds = find(strncmp('+', hit_data,1));
hit_data(inds) = num2cell(cellfun(@length, hit_data(inds)));

inds = find(strncmp('-', hit_data,1));
hit_data(inds) = num2cell(-cellfun(@length, hit_data(inds)));

inds = find(strcmp('a', hit_collection));
hit_strains_hap =  hit_strains(inds);
hit_data_hap = cell2mat(hit_data(inds));

inds = find(strcmp('a/a', hit_collection));
hit_strains_het = hit_strains(inds);
hit_data_het = cell2mat(hit_data(inds));

% If the same strain is present more than once, average its values
[hit_strains_hap, hit_data_hap] = grpstats(hit_data_hap, hit_strains_hap, {'gname','mean'});
[hit_strains_het, hit_data_het] = grpstats(hit_data_het, hit_strains_het, {'gname','mean'});

hit_strains = unique([hit_strains_hap; hit_strains_het]);
hit_data = nan(length(hit_strains),2);
[~,ind1,ind2] = intersect(hit_strains, hit_strains_hap);
hit_data(ind1,1) = hit_data_hap(ind2);
[~,ind1,ind2] = intersect(hit_strains, hit_strains_het);
hit_data(ind1,2) = hit_data_het(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [73; 414];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
zewaii_huang_2003.orfs = hit_strains;
zewaii_huang_2003.ph = hit_data_names;
zewaii_huang_2003.data = hit_data;
zewaii_huang_2003.dataset_ids = hit_data_ids;

%% Save

save('./zewaii_huang_2003.mat','zewaii_huang_2003');

%% Print out

fid = fopen('./zewaii_huang_2003.txt','w');
write_matrix_file(fid, zewaii_huang_2003.orfs, zewaii_huang_2003.ph, zewaii_huang_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(zewaii_huang_2003)
end

end

