%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
singh_babak_cowen_2012.pmid = 22751784;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(singh_babak_cowen_2012.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Singh-Babak2012_HIPHOPData.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(2:end,1);

hit_strains2 = regexpi(hit_strains, '[A-Z0-9\-]*(?=::)','match');
hit_strains2 = vertcat(hit_strains2{:});

% Get the data itself
hit_data = cell2mat(data(2:end,3:2:9));
hit_data_type = data(2:end, 17);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Normalize by DMSO
hit_data2 = hit_data - repmat(hit_data(:,1), 1, 4);
hit_data2(:,1) = [];

% Split het and hom
hit_data3 = nan(length(hit_strains), 6);

inds_hom = find(strcmp('hom', hit_data_type));
hit_data3(inds_hom,1:3) = hit_data2(inds_hom,:);
inds_het = find(strcmp('het', hit_data_type));
hit_data3(inds_het,4:6) = hit_data2(inds_het,:);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data3, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [5346 5348 5350 5343 5347 5349]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
singh_babak_cowen_2012.orfs = hit_strains;
singh_babak_cowen_2012.ph = hit_data_names;
singh_babak_cowen_2012.data = hit_data;
singh_babak_cowen_2012.dataset_ids = hit_data_ids;

%% Save

save('./singh_babak_cowen_2012.mat','singh_babak_cowen_2012');

%% Print out

fid = fopen('./singh_babak_cowen_2012.txt','w');
write_matrix_file(fid, singh_babak_cowen_2012.orfs, singh_babak_cowen_2012.ph, singh_babak_cowen_2012.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(singh_babak_cowen_2012)
end

end

