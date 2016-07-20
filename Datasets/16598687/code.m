%% van Voorst~Bradt, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
van_voorst_bradt_2006.pmid = 16598687;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(van_voorst_bradt_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES, hits] = read_data('textscan', './raw_data/hits_genenames.txt', '%s');

hits_orfs = translate(hits);

hits_scores = -ones(length(hits_orfs),1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [174];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

van_voorst_bradt_2006.orfs = hits_orfs;
van_voorst_bradt_2006.data = hits_scores;
van_voorst_bradt_2006.ph = hit_data_names;
van_voorst_bradt_2006.dataset_ids = hit_data_ids;

%% Save

save('./van_voorst_bradt_2006.mat','van_voorst_bradt_2006');

%% Print out

fid = fopen('./van_voorst_bradt_2006.txt','w');
write_matrix_file(fid, van_voorst_bradt_2006.orfs, van_voorst_bradt_2006.ph, van_voorst_bradt_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(van_voorst_bradt_2006)
end

end
