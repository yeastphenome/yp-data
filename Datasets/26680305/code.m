%% Wu~Hopper, 2015
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
wu_hopper_2015.pmid = 26680305;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(wu_hopper_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

%% First dataset
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS1.txt', '%s', 'delimiter', '\n');

% Split table based on spaces
C = cellfun(@strsplit, data, 'UniformOutput', 0);
data = C(cellfun(@length, C)>3);

% Get all the first columns
hit_strains1 = cellfun(@(v)v(1),data);

% Eliminate all white spaces & capitalize
hit_strains1 = clean_genename(hit_strains1);

% If in gene name form, transform into ORF name
[hit_strains1, translated, ambiguous] = translate(hit_strains1);

% If possible, fix the problem (typos, omissions etc.)
hit_strains1(ismember(hit_strains1, {'HXT12'})) = {'YIL170W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
hit_strains1(inds) = [];
data(inds) = [];

% Get indices where phenotype is by finding type:
ind = cellfun(@(s) find(~cellfun(@isempty, strfind(s, 'Type'))), data, 'Uni', 1);
phenotype = cellfun(@(c,idx)c(idx+1),data,num2cell(ind));
extra = find(contains(cellfun(@(c,idx)c(idx+2),data,num2cell(ind)), 'and'));

% Make the hit_data
hit_data1 = zeros(length(data), 6);
hit_data1(extra, 2) = 1; % for 1/2 + extra band
ind = find(~cellfun(@isempty, strfind(phenotype, '1/2'))); 
hit_data1(ind, 1) = 1; % for 1/2
ind = find(~cellfun(@isempty, strfind(phenotype, '2*'))); 
hit_data1(ind, 3) = 1; % for 2*
ind = find(~cellfun(@isempty, strfind(phenotype, '4/5/6'))); 
hit_data1(ind, 4) = 1; % for 4/5/6
ind = find(~cellfun(@isempty, strfind(phenotype, '7'))); 
hit_data1(ind, 5) = 1; % for 7
ind = find(~cellfun(@isempty, strfind(phenotype, '8'))); 
hit_data1(ind, 5) = 1; % for 8

%% Second dataset
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS2.txt', '%s', 'delimiter', '\n');

% Split table based on spaces
C = cellfun(@strsplit, data, 'UniformOutput', 0);
data = C(cellfun(@length, C)>3);

% Get all the first columns
hit_strains2 = cellfun(@(v)v(2),data);

% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
[hit_strains2, translated, ambiguous] = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
hit_strains2(inds) = [];
data(inds) = [];

% Get indices where phenotype is by finding type:
ind = cellfun(@(s) find(~cellfun(@isempty, strfind(s, 'Type'))), data, 'Uni', 1);
phenotype = cellfun(@(c,idx)c(idx+1),data,num2cell(ind));

% Make the hit_data
hit_data2 = zeros(length(data), 6);
ind = find(~cellfun(@isempty, strfind(phenotype, '1/2'))); 
hit_data2(ind, 1) = 1; % for 1/2
ind = find(~cellfun(@isempty, strfind(phenotype, '2*'))); 
hit_data2(ind, 3) = 1; % for 2*
ind = find(~cellfun(@isempty, strfind(phenotype, '4/5/6'))); 
hit_data2(ind, 4) = 1; % for 4/5/6
ind = find(~cellfun(@isempty, strfind(phenotype, '7'))); 
hit_data2(ind, 5) = 1; % for 7
ind = find(~cellfun(@isempty, strfind(phenotype, '8'))); 
hit_data2(ind, 5) = 1; % for 8

%% Second dataset
[FILENAMES{end+1}, data] = read_data('textread','./raw_data/TableS3.txt', '%s', 'delimiter', '\n');

% Split table based on spaces
C = cellfun(@strsplit, data, 'UniformOutput', 0);
data = C(cellfun(@length, C)>3);

% Get all the first columns
hit_strains3 = cellfun(@(v)v(2),data);

% Eliminate all white spaces & capitalize
hit_strains3 = clean_orf(hit_strains3);

% If in gene name form, transform into ORF name
[hit_strains3, translated, ambiguous] = translate(hit_strains3);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains3));
hit_strains3(inds) = [];
data(inds) = [];

% Get indices where phenotype is by finding type:
ind = cellfun(@(s) find(~cellfun(@isempty, strfind(s, 'Type'))), data, 'Uni', 1);
phenotype = cellfun(@(c,idx)c(idx+1),data,num2cell(ind));

% Make the hit_data
hit_data3 = zeros(length(data), 6);
ind = find(~cellfun(@isempty, strfind(phenotype, '1/2'))); 
hit_data3(ind, 1) = 1; % for 1/2
ind = find(~cellfun(@isempty, strfind(phenotype, '2*'))); 
hit_data3(ind, 3) = 1; % for 2*
ind = find(~cellfun(@isempty, strfind(phenotype, '4/5/6'))); 
hit_data3(ind, 4) = 1; % for 4/5/6
ind = find(~cellfun(@isempty, strfind(phenotype, '7'))); 
hit_data3(ind, 5) = 1; % for 7
ind = find(~cellfun(@isempty, strfind(phenotype, '8'))); 
hit_data3(ind, 5) = 1; % for 8

%% Combine hit datasets
hit_data = [hit_data1; hit_data2; hit_data3];
hit_strains = [hit_strains1; hit_strains2; hit_strains3];
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4952; 11779; 11780; 11781; 11782; 11783];

%% Tested strains 

% Load tested strains
[FILENAMES{end+1}, tested_strains1] = read_data('xlsread','./raw_data/Boone ts collection Vs. Hieter ts collection.xlsx', 'Boone only');
[FILENAMES{end+1}, tested_strains2] = read_data('xlsread','./raw_data/Boone ts collection Vs. Hieter ts collection.xlsx', 'Hieter only');

% Make one file
tested_strains = [tested_strains1; tested_strains2];

% Eliminate all white spaces & capitalize
tested_strains = clean_orf(tested_strains);

% If in gene name form, transform into ORF name
[tested_strains, ~, ~] = translate(tested_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
tested_strains(inds) = [];

% Finally, take the unique set
tested_strains = unique(tested_strains);

% Make sure the that all the hits are part of the tested set
[missing,~] = setdiff(hit_strains, tested_strains); % l09 missing

% If it seems reasonable, add the missing hits to the list of tested strains
tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
wu_hopper_2015.orfs = tested_strains;
wu_hopper_2015.ph = hit_data_names;
wu_hopper_2015.data = zeros(length(wu_hopper_2015.orfs),length(wu_hopper_2015.ph));
wu_hopper_2015.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, wu_hopper_2015.orfs);
wu_hopper_2015.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./wu_hopper_2015.mat','wu_hopper_2015');

%% Print out

fid = fopen('./wu_hopper_2015.txt','w');
write_matrix_file(fid, wu_hopper_2015.orfs, wu_hopper_2015.ph, wu_hopper_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(wu_hopper_2015)
end

end