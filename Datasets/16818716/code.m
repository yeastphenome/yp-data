%% Lam~Conibear, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
% TESTED = not available
lam_conibear_2006.pmid = 16818716;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(lam_conibear_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load data
[FILENAMES, C] = read_data('textscan', './raw_data/hits_genes_data.txt', '%s\t%.3f\n');

hits_genes = C{1};
hits_data = C{2};

hits_genes = clean_genename(hits_genes);

hits_orfs = translate(hits_genes);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [138];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
lam_conibear_2006.orfs = hits_orfs;
lam_conibear_2006.ph = hit_data_names;
lam_conibear_2006.data = hits_data;
lam_conibear_2006.dataset_ids = hit_data_ids;

%% Save

save('./lam_conibear_2006.mat','lam_conibear_2006');

%% Print out

fid = fopen('./lam_conibear_2006.txt','w');
write_matrix_file(fid, lam_conibear_2006.orfs, lam_conibear_2006.ph, lam_conibear_2006.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(lam_conibear_2006)
end

end
