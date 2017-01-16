%% Gustavsson~Ronne, 2008
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
gustavsson_ronne_2008.pmid = 18612650;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(gustavsson_ronne_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table1.xlsx');

hit_strains = {};
hit_data = {};

% Get the list of ORFs and the correponding data 
for i = 1:length(data)
    if isnan(data{i})
        continue
    end
    C = strsplit(data{i}, ' ');
    ind1 = regexp(C, '\d$');
    index = (cellfun(@(x) x==1, ind1, 'UniformOutput', 0));
    indx = find(cellfun(@any, index));
    if isempty(indx)
        continue
    end
    hit_strains = [hit_strains; C{1}];
    hit_data = [hit_data; C{indx(1)}];
end

hit_data = cell2mat(hit_data);

% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [105];

%% Prepare final dataset
% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
gustavsson_ronne_2008.orfs = hit_strains;
gustavsson_ronne_2008.ph = hit_data_names;
gustavsson_ronne_2008.data = hit_data;
gustavsson_ronne_2008.dataset_ids = hit_data_ids;

%% Save

save('./gustavsson_ronne_2008.mat','gustavsson_ronne_2008');

%% Print out

fid = fopen('./gustavsson_ronne_2008.txt','w');
write_matrix_file(fid, gustavsson_ronne_2008.orfs, gustavsson_ronne_2008.ph, gustavsson_ronne_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(gustavsson_ronne_2008)
end

end

