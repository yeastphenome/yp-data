%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
huang_kowalski_2005.pmid = 15994967;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(huang_kowalski_2005.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/hit_list.txt','ReadVariableNames', false);

% Get the list of ORFs and the correponding data 
hit_strains = data.Var1;

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains = unique(hit_strains);
hit_data = ones(size(hit_strains));


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [183];


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
huang_kowalski_2005.orfs = hit_strains;
huang_kowalski_2005.ph = hit_data_names;
huang_kowalski_2005.data = hit_data;
huang_kowalski_2005.dataset_ids = hit_data_ids;


%% Save

save('./huang_kowalski_2005.mat','huang_kowalski_2005');

%% Print out

fid = fopen('./huang_kowalski_2005.txt','w');
write_matrix_file(fid, huang_kowalski_2005.orfs, huang_kowalski_2005.ph, huang_kowalski_2005.data);
fclose(fid);

end

