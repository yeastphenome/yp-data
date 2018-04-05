%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

% Data obtained from R by running:
% > load(file="sup_file_Met.rda")
% > write.table(tfit.summary, file="Met_tfit_summary.txt", sep='\t', quote=FALSE, row.names=FALSE)

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
salignon_yvert_2018.pmid = 29507053;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(salignon_yvert_2018.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data (part 1)

[FILENAMES{end+1}, data1] = read_data('read_matrix_file','./raw_data/Salt_tfit_summary.txt', 2, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1.labels_row(:,1);

% Get the data itself
hit_data1 = data1.data(:,[1,2,6,7,8,9,10]);
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

% If possible, fix the problem (typos, omissions etc.)
for i = 1 : length(inds)
    hit_strains1{inds(i)} = hit_strains1{inds(i)}(1:end-1);
end

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids1 = [16167; 16168; 16169; 16175; 16176; 16177; 16178];


%% Load the data (part 2)

[FILENAMES{end+1}, data2] = read_data('read_matrix_file','./raw_data/Met_tfit_summary.txt', 2, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2.labels_row(:,1);

% Get the data itself
hit_data2 = data2.data(:,[1,2,6,7,8,9,10]);
   
% Eliminate all white spaces & capitalize
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If possible, fix the problem (typos, omissions etc.)
for i = 1 : length(inds)
    hit_strains2{inds(i)} = hit_strains2{inds(i)}(1:end-1);
end

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [16170; 16171; 16172; 16179; 16180; 16181; 16182];

%% Prepare final dataset

hit_data_ids = [hit_data_ids1; hit_data_ids2];

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains), length(hit_data_ids));

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1:7) = hit_data1(ind1,:);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,8:end) = hit_data2(ind1,:);

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
salignon_yvert_2018.orfs = hit_strains;
salignon_yvert_2018.ph = hit_data_names;
salignon_yvert_2018.data = hit_data;
salignon_yvert_2018.dataset_ids = hit_data_ids;

%% Save

save('./salignon_yvert_2018.mat','salignon_yvert_2018');

%% Print out

fid = fopen('./salignon_yvert_2018.txt','w');
write_matrix_file(fid, salignon_yvert_2018.orfs, salignon_yvert_2018.ph, salignon_yvert_2018.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(salignon_yvert_2018)
end

end

