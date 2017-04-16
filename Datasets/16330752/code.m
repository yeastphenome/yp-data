%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
proszynski_walch_2005.pmid = 16330752;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(proszynski_walch_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/09107Table2.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,2);
phenotypes = data(:,1);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
phenotypes(inds) = [];

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds)); 

hit_strains(inds) = [];
phenotypes(inds) = [];

inds1 = find(ismember(phenotypes, {'I','I and II','II and I'}));
inds2 = find(ismember(phenotypes, {'II','I and II','II and I'}));

hit_data = zeros(length(hit_strains),2);
hit_data(inds1,1) = 1;
hit_data(inds2,2) = 1;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [180; 5657];


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
proszynski_walch_2005.orfs = hit_strains;
proszynski_walch_2005.ph = hit_data_names;
proszynski_walch_2005.data = hit_data;
proszynski_walch_2005.dataset_ids = hit_data_ids;

%% Save

save('./proszynski_walch_2005.mat','proszynski_walch_2005');

%% Print out

fid = fopen('./proszynski_walch_2005.txt','w');
write_matrix_file(fid, proszynski_walch_2005.orfs, proszynski_walch_2005.ph, proszynski_walch_2005.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(proszynski_walch_2005)
end

end

