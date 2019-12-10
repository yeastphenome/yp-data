%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
albulescu_pleiss_2012.pmid = 22479188;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(albulescu_pleiss_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/hits.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

% Get the data itself
hit_data = data(2:end,3:7);
hit_data = cell2mat(hit_data);
hit_data(isnan(hit_data)) = 1;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16311; 16315; 16312; 16314; 16313];

%% Tested strains (only if the dataset is not quantitative and the tested strains are provided separately)

% Load tested strains
[FILENAMES{end+1}, tested_strains] = read_data('textread','./raw_data/no_hits.txt','%s');

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
tested_strains = translate(tested_strains);

% If possible, fix the typo
tested_strains(ismember(tested_strains, {'YLR287-A'})) = {'YLR287C-A'};
tested_strains(ismember(tested_strains, {'YCRO54C'})) = {'YCR054C'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains);
disp(missing);


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
albulescu_pleiss_2012.orfs = tested_strains;
albulescu_pleiss_2012.ph = hit_data_names;
albulescu_pleiss_2012.data = ones(length(albulescu_pleiss_2012.orfs),length(albulescu_pleiss_2012.ph));
albulescu_pleiss_2012.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, albulescu_pleiss_2012.orfs);
albulescu_pleiss_2012.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./albulescu_pleiss_2012.mat','albulescu_pleiss_2012');

%% Print out

fid = fopen('./albulescu_pleiss_2012.txt','w');
write_matrix_file(fid, albulescu_pleiss_2012.orfs, albulescu_pleiss_2012.ph, albulescu_pleiss_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(albulescu_pleiss_2012)
end

end

