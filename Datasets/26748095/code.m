%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
simoneau_wurtele_2016.pmid = 26748095;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(simoneau_wurtele_2016.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/nar-02322-x-2015-File011.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

t = split(hit_strains, '::');
hit_strains = t(:,1);

% Get the data itself (taking the opposite because the Z-scores were
% calculated as untreated/treated)
hit_data = -cell2mat(data(3:end,3:6));
   
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
hit_data_ids = [16215; 16214; 16215; 16214];
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
simoneau_wurtele_2016.orfs = hit_strains;
simoneau_wurtele_2016.ph = hit_data_names;
simoneau_wurtele_2016.data = hit_data;
simoneau_wurtele_2016.dataset_ids = hit_data_ids;

%% Save

save('./simoneau_wurtele_2016.mat','simoneau_wurtele_2016');

%% Print out

fid = fopen('./simoneau_wurtele_2016.txt','w');
write_matrix_file(fid, simoneau_wurtele_2016.orfs, simoneau_wurtele_2016.ph, simoneau_wurtele_2016.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(simoneau_wurtele_2016)
end

end

