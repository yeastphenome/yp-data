%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
willingham_muchowski_2003.pmid = 14657499;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(willingham_muchowski_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/alpha-sinuclein.txt', 'delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data.Strain;
tmp = regexp(hit_strains1, ' ', 'split');
hit_strains1 = cellfun(@(c) c{2}, tmp, 'UniformOutput', false);

% Get the data itself
hit_data1 = -ones(size(hit_strains1)); % if the dataset is binary
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_genename(hit_strains1);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids1 = [202];

%% Load the data 2

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/huntingtin.txt', 'delimiter','\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data.Strain;
tmp = regexp(hit_strains2, ' ', 'split');
hit_strains2 = cellfun(@(c) c{2}, tmp, 'UniformOutput', false);

% Get the data itself
hit_data2 = -ones(size(hit_strains2)); % if the dataset is binary
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_genename(hit_strains2);

hit_strains2(strcmp('pc16', hit_strains2)) = {'pcl6'};
hit_strains2(strcmp('mrp11', hit_strains2)) = {'mrpl1'};

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [70];

%% Merge the 2 datasets

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = zeros(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1) = -1;
[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,2) = -1;

hit_data_ids = [hit_data_ids1; hit_data_ids2];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
willingham_muchowski_2003.orfs = hit_strains;
willingham_muchowski_2003.ph = hit_data_names;
willingham_muchowski_2003.data = hit_data;
willingham_muchowski_2003.dataset_ids = hit_data_ids;

%% Save

save('./willingham_muchowski_2003.mat','willingham_muchowski_2003');

%% Print out

fid = fopen('./willingham_muchowski_2003.txt','w');
write_matrix_file(fid, willingham_muchowski_2003.orfs, willingham_muchowski_2003.ph, willingham_muchowski_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(willingham_muchowski_2003)
end

end

