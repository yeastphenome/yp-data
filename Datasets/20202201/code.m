%% Batova~Schuller, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
batova_schuller_2010.pmid = 20202201;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(batova_schuller_2010.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/1471-2164-11-153-S1.xlsx');

% Clean through the data
data = data(21:end, 7);
data(find(cellfun(@isempty, data))) = [];
indx = find(cell2mat(arrayfun(@is_orf, data, 'UniformOutput', false)));

% Isolate the hit strains
hit_strains = data(indx);

% Clean hit_strains
hit_strains = clean_orf(hit_strains);

% Make a data matrix
hit_data = zeros(length(hit_strains), 2);

% Fill in the data matrix
for i = 1:length(indx)
    if ~cellfun(@isempty, strfind(data(indx(i)+2), '.'))
        hit_data(i, 1) = -2;
    elseif ~cellfun(@isempty, strfind(data(indx(i)+2), '+'))
        hit_data(i, 1) = 1;
    elseif ~cellfun(@isempty, strfind(data(indx(i)+2), 'sl'))
        hit_data(i, 1) = -1; 
    end
    
    if ~cellfun(@isempty, strfind(data(indx(i)+3), '.'))
        hit_data(i, 2) = -2;
    elseif ~cellfun(@isempty, strfind(data(indx(i)+3), '+'))
        hit_data(i, 2) = 1;
    elseif ~cellfun(@isempty, strfind(data(indx(i)+3), 'sl'))
        hit_data(i, 2) = -1; 
    end
end

% Clean up orf names and remove repeats
hit_strains = clean_orf(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [153; 437];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
batova_schuller_2010.orfs = hit_strains;
batova_schuller_2010.ph = hit_data_names;
batova_schuller_2010.data = hit_data;
batova_schuller_2010.dataset_ids = hit_data_ids;

%% Save

save('./batova_schuller_2010.mat','batova_schuller_2010');

%% Print out

fid = fopen('./batova_schuller_2010.txt','w');
write_matrix_file(fid, batova_schuller_2010.orfs, batova_schuller_2010.ph, batova_schuller_2010.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(batova_schuller_2010)
end

end

