%% Serrano~Arino, 2004
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
serrano_arino_2004.pmid = 14993228;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(serrano_arino_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/serrano_arino_2004.xlsx', 'Sheet1');

% Eliminate white spaces before/after ORF
hit_genenames = clean_genename(data.raw(:,1));

% Translate genenames to ORF
hit_strains = translate(hit_genenames);

hit_data = cell2mat(data.raw(:,2));

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

% Flip the scores such that -5 is the most sensitive and -1 is the least
% sensitive
hit_data = hit_data - 6;


% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [187];

%% Load tested genes
[FILENAMES{end+1}, tested.raw] = read_data('xlsread','./raw_data/BY4741.xlsx', 'Tabelle1');
tested_strains = clean_orf(tested.raw(2:end,2));

% Find anything that doesn't look like an ORF
inds = find(~is_orf(tested_strains));
disp(tested_strains(inds));  

tested_strains = unique(tested_strains);

% Check if all the hits are in the tested space
[missing,inds] = setdiff(hit_strains, tested_strains);
tested_strains = [tested_strains; missing];   % addding 3 mising ORFs

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
serrano_arino_2004.orfs = tested_strains;
serrano_arino_2004.ph = hit_data_names;
serrano_arino_2004.data = zeros(length(serrano_arino_2004.orfs),length(serrano_arino_2004.ph));
serrano_arino_2004.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hit_strains, serrano_arino_2004.orfs);
serrano_arino_2004.data(ind2,:) = hit_data(ind1,:);

%% Save

save('./serrano_arino_2004.mat','serrano_arino_2004');

%% Print out

fid = fopen('./serrano_arino_2004.txt','w');
write_matrix_file(fid, serrano_arino_2004.orfs, serrano_arino_2004.ph, serrano_arino_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(serrano_arino_2004)
end

end

