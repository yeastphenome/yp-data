%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
markovich_osherov_2004.pmid = 15388447;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(markovich_osherov_2004.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('textscan','./raw_data/hits.txt','%s','delimiter',',');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1);
ind1 = find(strcmp('Sensitive', hit_strains));
ind2 = find(strcmp('Resistant', hit_strains));

hit_strains1 = hit_strains(ind1+1:ind2-1);
hit_strains2 = hit_strains(ind2+1:end);

hit_strains = [hit_strains1; hit_strains2];
hit_data = [-ones(length(hit_strains1),1); ones(length(hit_strains2),1)];

hit_strains(strcmp('SMI1/KNR4', hit_strains)) = {'SMI1'};
hit_strains(strcmp('CHS4/SKT5', hit_strains)) = {'CHS4'};
   
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
hit_data_ids = [189];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
markovich_osherov_2004.orfs = hit_strains;
markovich_osherov_2004.ph = hit_data_names;
markovich_osherov_2004.data = hit_data;
markovich_osherov_2004.dataset_ids = hit_data_ids;

%% Save

save('./markovich_osherov_2004.mat','markovich_osherov_2004');

%% Print out

fid = fopen('./markovich_osherov_2004.txt','w');
write_matrix_file(fid, markovich_osherov_2004.orfs, markovich_osherov_2004.ph, markovich_osherov_2004.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(markovich_osherov_2004)
end

end

