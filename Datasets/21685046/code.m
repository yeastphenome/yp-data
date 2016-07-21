%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yu_vitek_2011.pmid = 21685046;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yu_vitek_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data (hap)

[FILENAMES{end+1}, data] = read_data('read_matrix_file','./raw_data/PiiMS_Dataset_Full_genome_knockout_haploid_1469056539143/Dataset_387_Aggregate_1299870873410.csv', 'delimiter', ',');

% Get the list of ORFs and the correponding data 
hit_strains = data.labels_row;

% Get the data itself
hit_data = data.data;
   
% Eliminate all white spaces & capitalize
t = cellfun(@(x) regexp(x, '_', 'split'), hit_strains,'UniformOutput',0);
hit_strains = cellfun(@(x) x{1}, t, 'UniformOutput', 0);
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
hit_strains(strcmp('YLR287-A', hit_strains)) = {'YLR287C-A'};
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [4797 4798 4799 4800 4801 4802 4803 4804 4805 4806 4807 4808 4809 4810]';

%% Load the data (het)

[FILENAMES{end+1}, data2] = read_data('read_matrix_file','./raw_data/PiiMS_Dataset_Full_genome_knockout_diploid_1469056540949/Dataset_390_Aggregate_1299877732527.csv', 'delimiter', ',');

% Get the list of ORFs and the correponding data 
hit_strains2 = data2.labels_row;

% Get the data itself
hit_data2 = data2.data;
   
% Eliminate all white spaces & capitalize
t = cellfun(@(x) regexp(x, '_', 'split'), hit_strains2,'UniformOutput',0);
hit_strains2 = cellfun(@(x) x{1}, t, 'UniformOutput', 0);
hit_strains2 = clean_orf(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [4811:4824]';

%%

hit_strains3 = unique([hit_strains; hit_strains2]);
hit_data3 = nan(length(hit_strains3),28);
[~,ind1,ind2] = intersect(hit_strains, hit_strains3);
hit_data3(ind2,1:14) = hit_data(ind1,:);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains3);
hit_data3(ind2,15:28) = hit_data2(ind1,:);

hit_data_ids3 = [hit_data_ids; hit_data_ids2];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids3);
hit_data_names3 = cell(size(hit_data_ids3));
hit_data_names3(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yu_vitek_2011.orfs = hit_strains3;
yu_vitek_2011.ph = hit_data_names3;
yu_vitek_2011.data = hit_data3;
yu_vitek_2011.dataset_ids = hit_data_ids3;

%% Save

save('./yu_vitek_2011.mat','yu_vitek_2011');

%% Print out

fid = fopen('./yu_vitek_2011.txt','w');
write_matrix_file(fid, yu_vitek_2011.orfs, yu_vitek_2011.ph, yu_vitek_2011.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yu_vitek_2011)
end

end

