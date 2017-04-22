%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
giaever_davis_2004.pmid = 14718668;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(giaever_davis_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the Hom data

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/analyzed_data_hom.txt', 2, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hom_strains = data.labels_row(:,1);

% Get the data itself
hom_data = data.data;
   
% Eliminate all white spaces & capitalize
hom_strains = clean_orf(hom_strains);

% If in gene name form, transform into ORF name
hom_strains = translate(hom_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hom_strains));
disp(hom_strains(inds));  

% If the same strain is present more than once, average its values
[hom_strains, hom_data] = grpstats(hom_data, hom_strains, {'gname','mean'});

% Clean-up
[FILENAMES{end+1}, eliminated1] = read_data('textscan','./raw_data/background_tags_base_hom.txt', '%s');

eliminated1 = clean_orf(eliminated1);
eliminated1 = translate(eliminated1);

[~,~,ind2] = intersect(eliminated1, hom_strains);
hom_strains(ind2) = [];
hom_data(ind2,:) = [];

% Remove all ORFs that only have 0s because >90% of them are essential
% genes
nZ = sum(hom_data == 0,2);
inds = find(nZ == size(hom_data,2));
hom_strains(inds) = [];
hom_data(inds,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hom_data_ids = [5316 5316 5316 5317 5317 5318 5318 5313 5313]';
[hom_data_ids, hom_data] = grpstats(hom_data', hom_data_ids, {'gname','mean'});
hom_data = hom_data';
hom_data_ids = cellfun(@str2num, hom_data_ids);

%% Load Het data

[FILENAMES{end+1}, data2] = read_data('read_matrix_file','./raw_data/analyzed_data_het.txt', 2, 1);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
het_strains = data2.labels_row(:,1);

% Get the data itself
het_data = data2.data;
   
% Eliminate all white spaces & capitalize
het_strains = clean_orf(het_strains);

% If in gene name form, transform into ORF name
het_strains = translate(het_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(het_strains));
disp(het_strains(inds));  

% If the same strain is present more than once, average its values
[het_strains, het_data] = grpstats(het_data, het_strains, {'gname','mean'});

% Clean-up
[FILENAMES{end+1}, eliminated_het1] = read_data('textscan','./raw_data/background_tags_base.txt', '%s');

eliminated_het1 = clean_orf(eliminated_het1);
eliminated_het1 = translate(eliminated_het1);

[~,~,ind2] = intersect(eliminated_het1, het_strains);
het_strains(ind2) = [];
het_data(ind2,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
het_data_ids = [5312 5312 5312 5312 5312 5314 5314 5310 5315 5315 5311 5337 5337 5337 5337 5323 5323 5322 5322 5322]';
het_data_ids = [het_data_ids; 5321; 5321; 5321; 5319; 5319; 5320; 5320; 5338; 5339; 5339; 5339; 5339];
het_data_ids = [het_data_ids; 5334; 5335; 5335; 5336; 5336; 5333; 5333; 5330; 5331; 5331; 5332; 5332; 5325; 5325; 5325; 5324];
het_data_ids = [het_data_ids; 484; 484; 484; 5308; 5308; 5308; 5308; 5308; 5308; 5308; 5308; 5308; 5309];
het_data_ids = [het_data_ids; 5326; 5326; 5327; 5327; 5327; 5327; 5328; 5328; 5329; 5329];
[het_data_ids, het_data] = grpstats(het_data', het_data_ids, {'gname','mean'});
het_data = het_data';
het_data_ids = cellfun(@str2num, het_data_ids);

%% Prepare final dataset

hit_strains = unique([hom_strains; het_strains]);
hit_data_ids = [hom_data_ids; het_data_ids];
hit_data = nan(length(hit_strains), length(hit_data_ids));

[~,ind1,ind2] = intersect(hit_strains, hom_strains);
hit_data(ind1,1:length(hom_data_ids)) = hom_data(ind2,:);

[~,ind1,ind2] = intersect(hit_strains, het_strains);
hit_data(ind1,length(hom_data_ids)+1:end) = het_data(ind2,:);

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
giaever_davis_2004.orfs = hit_strains;
giaever_davis_2004.ph = hit_data_names;
giaever_davis_2004.data = -hit_data;  % taking the opposite to make sure that low number indicate low growth and viceversa
giaever_davis_2004.dataset_ids = hit_data_ids;

%% Save

save('./giaever_davis_2004.mat','giaever_davis_2004');

%% Print out

fid = fopen('./giaever_davis_2004.txt','w');
write_matrix_file(fid, giaever_davis_2004.orfs, giaever_davis_2004.ph, giaever_davis_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(giaever_davis_2004)
end

end

