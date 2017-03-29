%% Alto~Dixon, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
alto_dixon_2006.pmid = 16413487;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(alto_dixon_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data

[FILENAMES{end+1}, hits_gn] = read_data('textread','./raw_data/hits.txt', '%s');
hits_data = ones(size(hits_gn));

hits_gn = clean_genename(hits_gn);

hits_orfs = translate(hits_gn);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [561];

%% Load tested

[FILENAMES{end+1}, tested_orfs] = read_data('textread','./raw_data/FG_array_genes.txt', '%s');

tested_orfs = clean_orf(tested_orfs);

[missing, ix] = setdiff(hits_orfs, tested_orfs);

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
alto_dixon_2006.orfs = tested_orfs;
alto_dixon_2006.ph = hit_data_names;
alto_dixon_2006.data = zeros(length(alto_dixon_2006.orfs),length(alto_dixon_2006.ph));
alto_dixon_2006.dataset_ids = hit_data_ids;

[~,ind1,ind2] = intersect(hits_orfs, alto_dixon_2006.orfs);
alto_dixon_2006.data(ind2,:) = hits_data(ind1,:);

%% Save

save('./alto_dixon_2006.mat','alto_dixon_2006');

%% Print out

fid = fopen('./alto_dixon_2006.txt','w');
write_matrix_file(fid, alto_dixon_2006.orfs, alto_dixon_2006.ph, alto_dixon_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(alto_dixon_2006)
end

end
