%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
papic_rapaport_2013.pmid = 23959800;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(papic_rapaport_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textread','./raw_data/hits.txt', '%s %s %s','delimiter', '\t');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data{1};

% Get the data itself
hit_data = cat(2, data{2}, data{3});
   
% Eliminate all white spaces & capitalize
hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

hit_strains(inds) = [];
hit_data(inds,:) = [];

hit_data = cellfun(@str2num, hit_data);

% flipping the sign so that, when the "non-hit" strains are set to 0, 
% the overall values reflect the extent to which the phenotype (protein insertion into MOM)
% is expressed: low values = defect in phenotype (incorrect insertion), 
% high values = no phenotype (correct insertion)
hit_data = -hit_data; 

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

hit_data = nanmean(hit_data, 2);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11867];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
papic_rapaport_2013.orfs = hit_strains;
papic_rapaport_2013.ph = hit_data_names;
papic_rapaport_2013.data = hit_data;
papic_rapaport_2013.dataset_ids = hit_data_ids;

%% Save

save('./papic_rapaport_2013.mat','papic_rapaport_2013');

%% Print out

fid = fopen('./papic_rapaport_2013.txt','w');
write_matrix_file(fid, papic_rapaport_2013.orfs, papic_rapaport_2013.ph, papic_rapaport_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(papic_rapaport_2013)
end

end

