%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
ericson_nislow_2008.pmid = 18688276;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(ericson_nislow_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pgen.1000151.s003.xlsx', 'Sheet1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(4:end,1);
hit_hom_het = data(4:end,95);

hit_screens = data(3,3:88)';

% Get the data itself
hit_data = cell2mat(data(4:end,3:88));
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% Separate hom screens and het screens

hit_strains2 = hit_strains;
hit_screens2 = [hit_screens; hit_screens];
hit_data2 = nan(length(hit_strains2), length(hit_screens2));

inds = find(strcmp('hom', hit_hom_het));
hit_data2(inds,1:length(hit_screens)) = hit_data(inds,:);

inds = find(strcmp('het', hit_hom_het));
hit_data2(inds,length(hit_screens)+1:end) = hit_data(inds,:);

hit_strains = hit_strains2;
hit_screens = hit_screens2;
hit_data = hit_data2;

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
[FILENAMES{end+1}, data] = read_data('xlsread','./extras/drug_dataset.xlsx', 'Sheet1');
drugs = data(2:87,1);
hom_id = cell2mat(data(2:87,3));
het_id = cell2mat(data(2:87,4));

[~,inds] = ismember(hit_screens(1:length(hit_screens)/2), drugs);
hit_data_ids_hom = hom_id(inds);
hit_data_ids_het = het_id(inds);

hit_data_ids = [hit_data_ids_hom; hit_data_ids_het];
[hit_data_ids, hit_data] = grpstats(hit_data', hit_data_ids, {'gname','mean'});
hit_data = hit_data';
hit_data_ids = cellfun(@str2num, hit_data_ids);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
ericson_nislow_2008.orfs = hit_strains;
ericson_nislow_2008.ph = hit_data_names;
ericson_nislow_2008.data = hit_data;
ericson_nislow_2008.dataset_ids = hit_data_ids;

%% Save

save('./ericson_nislow_2008.mat','ericson_nislow_2008');

%% Print out

fid = fopen('./ericson_nislow_2008.txt','w');
write_matrix_file(fid, ericson_nislow_2008.orfs, ericson_nislow_2008.ph, ericson_nislow_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(ericson_nislow_2008)
end

end

