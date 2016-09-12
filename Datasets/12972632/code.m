%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
huang_kolodner_2003.pmid = 12972632;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(huang_kolodner_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/hits.txt', 'delimiter','\t','ReadVariableNames',false);

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.Var1;

% Get the data itself
hit_data = data.Var2;

inds = find(strcmp('Wild type', hit_strains));
hit_data = hit_data ./ hit_data(inds);

hit_strains(inds) = [];
hit_data(inds) = [];
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

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
hit_data_ids = [74];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
huang_kolodner_2003.orfs = hit_strains;
huang_kolodner_2003.ph = hit_data_names;
huang_kolodner_2003.data = hit_data;
huang_kolodner_2003.dataset_ids = hit_data_ids;

%% Save

save('./huang_kolodner_2003.mat','huang_kolodner_2003');

%% Print out

fid = fopen('./huang_kolodner_2003.txt','w');
write_matrix_file(fid, huang_kolodner_2003.orfs, huang_kolodner_2003.ph, huang_kolodner_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(huang_kolodner_2003)
end

end

