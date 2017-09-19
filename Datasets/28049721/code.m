%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
yofe_thoms_2017.pmid = 28049721;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(yofe_thoms_2017.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/TablesS1-S7.xlsx', 'Table S1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(3:end,1);

% Get the data itself
hit_data = data(3:end,[5 8]);
   
% Eliminate all white spaces & capitalize
hit_strains = clean_orf(hit_strains);

hit_strains(strcmp('YOLO57W', hit_strains)) = {'YOL057W'};
hit_strains(strcmp('YOLO62C', hit_strains)) = {'YOL062C'};
hit_strains(strcmp('YKLO72W', hit_strains)) = {'YKL072W'};
hit_strains(strcmp('YJL206-A', hit_strains)) = {'YJL206C-A'};
hit_strains(strcmp('YLR287-A', hit_strains)) = {'YLR287C-A'};

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));  

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11818 11819]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
yofe_thoms_2017.orfs = hit_strains;
yofe_thoms_2017.ph = hit_data_names;
yofe_thoms_2017.data = hit_data;
yofe_thoms_2017.dataset_ids = hit_data_ids;

% Data contains both deletions of nonessential genes and damp alleles of
% essential genes. Must remove the essential genes (no other way available)

load essential_genes_151215.mat
[~,ind1,ind2] = intersect(essential_genes, yofe_thoms_2017.orfs);
yofe_thoms_2017.orfs(ind2) = [];
yofe_thoms_2017.data(ind2,:) = [];


%% Save

save('./yofe_thoms_2017.mat','yofe_thoms_2017');

%% Print out

fid = fopen('./yofe_thoms_2017.txt','w');
write_matrix_file(fid, yofe_thoms_2017.orfs, yofe_thoms_2017.ph, yofe_thoms_2017.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(yofe_thoms_2017)
end

end

