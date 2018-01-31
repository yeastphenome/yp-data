%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
richie_hoepfner_2013.pmid = 23478965;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(richie_hoepfner_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the HOM data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/HOP-exp-scores-annotation.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data.SystematicNames;

% Get the data itself
hit_data = data.Score;

hit_compound_conc = data.Compound_conc;

% Make data into a matrix
hit_strains_u = unique(hit_strains);
hit_compound_conc_u = unique(hit_compound_conc);
hit_compound_conc_u(find(strcmp(hit_compound_conc_u, ''))) = [];
hit_data_matrix = nan(length(hit_strains_u), length(hit_compound_conc_u));

for i = 1 : length(hit_compound_conc_u)
    inds = find(strcmp(hit_compound_conc_u{i}, hit_compound_conc));
    if length(inds) > length(unique(hit_strains(inds)))
        i
    end
    [~,ind1,ind2] = intersect(hit_strains(inds), hit_strains_u);
    hit_data_matrix(ind2,i) = hit_data(inds(ind1));
end

% Eliminate all white spaces & capitalize
hit_strains_u = clean_orf(hit_strains_u);

% If in gene name form, transform into ORF name
hit_strains_u = translate(hit_strains_u);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_u));
disp(hit_strains_u(inds));  

hit_strains_u(inds) = [];
hit_data_matrix(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data_matrix, hit_strains_u, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16018, 16008, 16009, 16010, 16011, 16012, 16013, 16014, 16015, 16016, 16017]';

%% Load the HET data

[FILENAMES{end+1}, data2] = read_data('readtable','./raw_data/HIP-exp-scores-annotation.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains2 = data2.SystematicNames;

% Get the data itself
hit_data2 = data2.Score;

hit_compound_conc2 = data2.Compound_conc;

% Make data into a matrix
hit_strains_u2 = unique(hit_strains2);
hit_compound_conc_u2 = unique(hit_compound_conc2);
hit_compound_conc_u2(find(strcmp(hit_compound_conc_u2, ''))) = [];
hit_data_matrix2 = nan(length(hit_strains_u2), length(hit_compound_conc_u2));

for i = 1 : length(hit_compound_conc_u2)
    inds = find(strcmp(hit_compound_conc_u2{i}, hit_compound_conc2));
    if length(inds) > length(unique(hit_strains2(inds)))
        i
    end
    [~,ind1,ind2] = intersect(hit_strains2(inds), hit_strains_u2);
    hit_data_matrix2(ind2,i) = hit_data2(inds(ind1));
end

% Eliminate all white spaces & capitalize
hit_strains_u2 = clean_orf(hit_strains_u2);

% If in gene name form, transform into ORF name
hit_strains_u2 = translate(hit_strains_u2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains_u2));
disp(hit_strains_u2(inds));  

hit_strains_u2(inds) = [];
hit_data_matrix2(inds,:) = [];

% If the same strain is present more than once, average its values
[hit_strains2, hit_data2] = grpstats(hit_data_matrix2, hit_strains_u2, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids2 = [16029, 16019:16028]';

%% Prepare final dataset

% Merge HOM and HET

hit_strains_all = unique([hit_strains; hit_strains2]);
hit_data_ids_all = [hit_data_ids; hit_data_ids2];
hit_data_all = nan(length(hit_strains_all), length(hit_data_ids_all));

[~,ind1,ind2] = intersect(hit_strains, hit_strains_all);
hit_data_all(ind2,1:11) = hit_data(ind1,:);

[~,ind1,ind2] = intersect(hit_strains2, hit_strains_all);
hit_data_all(ind2,12:end) = hit_data2(ind1,:);

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids_all);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
richie_hoepfner_2013.orfs = hit_strains_all;
richie_hoepfner_2013.ph = hit_data_names;
richie_hoepfner_2013.data = hit_data_all;
richie_hoepfner_2013.dataset_ids = hit_data_ids_all;

%% Save

save('./richie_hoepfner_2013.mat','richie_hoepfner_2013');

%% Print out

fid = fopen('./richie_hoepfner_2013.txt','w');
write_matrix_file(fid, richie_hoepfner_2013.orfs, richie_hoepfner_2013.ph, richie_hoepfner_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(richie_hoepfner_2013)
end

end

