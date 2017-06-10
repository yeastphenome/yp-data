%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
choi_oshea_2017.pmid = 28520786;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(choi_oshea_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

filename = './raw_data/journal.pone.0176085.s004.xlsx';
[status, sheets] = xlsfinfo(filename);

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/journal.pone.0176085.s004.xlsx', sheets{1});
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/journal.pone.0176085.s004.xlsx', sheets{3});

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,1);
hit_strains2 = data2(2:end,1);

% Get the data itself
hit_data1 = data1(2:end,3); 
hit_data2 = data2(2:end,3); 
   
% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

hit_strains1(strcmp('YLR287-A', hit_strains1)) = {'YLR287C-A'};
hit_strains2(strcmp('YLR287-A', hit_strains2)) = {'YLR287C-A'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1) = hit_data1(ind1,1);

[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,2) = hit_data2(ind1,1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11807 11808]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
choi_oshea_2017.orfs = hit_strains;
choi_oshea_2017.ph = hit_data_names;
choi_oshea_2017.data = hit_data;
choi_oshea_2017.dataset_ids = hit_data_ids;

%% Save

save('./choi_oshea_2017.mat','choi_oshea_2017');

%% Print out

fid = fopen('./choi_oshea_2017.txt','w');
write_matrix_file(fid, choi_oshea_2017.orfs, choi_oshea_2017.ph, choi_oshea_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(choi_oshea_2017)
end

end

