%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
vizeacoumar_andrews_2010.pmid = 20065090;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(vizeacoumar_andrews_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/JCB_200909013_TableS3.xlsx','Single deletion screen');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

% Get the data itself
hit_data = cell2mat(data(3:end,2:91));
   
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
hit_data_ids = [16036:16125]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
vizeacoumar_andrews_2010.orfs = hit_strains;
vizeacoumar_andrews_2010.ph = hit_data_names;
vizeacoumar_andrews_2010.data = hit_data;
vizeacoumar_andrews_2010.dataset_ids = hit_data_ids;

%% Save

save('./vizeacoumar_andrews_2010.mat','vizeacoumar_andrews_2010');

%% Print out

fid = fopen('./vizeacoumar_andrews_2010.txt','w');
write_matrix_file(fid, vizeacoumar_andrews_2010.orfs, vizeacoumar_andrews_2010.ph, vizeacoumar_andrews_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(vizeacoumar_andrews_2010)
end

end

