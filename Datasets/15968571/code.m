%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
delia_hammond_2005.pmid = 15968571;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(delia_hammond_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('textscan','./raw_data/hits.txt', '%s %f');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1{1};

% Get the data itself
hit_data1 = data1{2};
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

%%

[FILENAMES{end+1}, data2] = read_data('textscan','./raw_data/hits2.txt', '%s %f');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2{1};

% Get the data itself
hit_data2 = data2{2};
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

%%

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = zeros(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1) = hit_data1(ind2);
[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,2) = hit_data2(ind2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [495; 5006];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
delia_hammond_2005.orfs = hit_strains;
delia_hammond_2005.ph = hit_data_names;
delia_hammond_2005.data = hit_data;
delia_hammond_2005.dataset_ids = hit_data_ids;

%% Save

save('./delia_hammond_2005.mat','delia_hammond_2005');

%% Print out

fid = fopen('./delia_hammond_2005.txt','w');
write_matrix_file(fid, delia_hammond_2005.orfs, delia_hammond_2005.ph, delia_hammond_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(delia_hammond_2005)
end

end

