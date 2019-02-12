%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
aiyar_steinmetz_2014.pmid = 25519239;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(aiyar_steinmetz_2014.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('readtable','./raw_data/ncomms6585-s3.txt');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
t = split(data.Strain,'::');
hit_strains = t(:,1);

% Get the data itself
hit_data = data.z_score;
zyg = data.Zygosity;
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_strains1 = hit_strains(strcmp('hom', zyg));
hit_data1 = hit_data(strcmp('hom', zyg));

hit_strains2 = hit_strains(strcmp('het', zyg));
hit_data2 = hit_data(strcmp('het', zyg));

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% Merge
hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = nan(length(hit_strains),2);

[~,ind1,ind2] = intersect(hit_strains, hit_strains1);
hit_data(ind1,1) = hit_data1(ind2);
[~,ind1,ind2] = intersect(hit_strains, hit_strains2);
hit_data(ind1,2) = hit_data2(ind2);


% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [16251; 16252];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
aiyar_steinmetz_2014.orfs = hit_strains;
aiyar_steinmetz_2014.ph = hit_data_names;
aiyar_steinmetz_2014.data = -hit_data;  % flip the sign since z-scores were calculated as (u_ctrl-u_treat)/s_ctrl
aiyar_steinmetz_2014.dataset_ids = hit_data_ids;

%% Save

save('./aiyar_steinmetz_2014.mat','aiyar_steinmetz_2014');

%% Print out

fid = fopen('./aiyar_steinmetz_2014.txt','w');
write_matrix_file(fid, aiyar_steinmetz_2014.orfs, aiyar_steinmetz_2014.ph, aiyar_steinmetz_2014.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(aiyar_steinmetz_2014)
end

end

