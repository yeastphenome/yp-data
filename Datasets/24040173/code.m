%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
novo_gonzalez_2013.pmid = 24040173;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(novo_gonzalez_2013.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

%% Phase 1
[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/Only original foldchanges/Direct comparison/Phase I/HOP_t0vsYPD10 (dir) PhI.xlsx', 'HOPt0vsYPD10_t0');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/Only original foldchanges/Direct comparison/Phase I/HOP_t0vsMS10 (dir) PhI.xlsx', 'HOPt0vsMS10_t0');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains1 = data1(2:end,1);
hit_data1 = data1(2:end,2);
essential1 = data1(2:end,5);

inds = find(strcmp('yes', essential1));
hit_strains1(inds) = [];
hit_data1(inds) = [];

hit_strains2 = data2(2:end,1);
hit_data2 = data2(2:end,2);
essential2 = data2(2:end,5);

inds = find(strcmp('yes', essential2));
hit_strains2(inds) = [];
hit_data2(inds) = [];
  
% Get ORFs
hit_strains1 = cellfun(@(x) regexp(x, '([\w\-])*(?=::)', 'match'), hit_strains1);
hit_strains2 = cellfun(@(x) regexp(x, '([\w\-])*(?=::)', 'match'), hit_strains2);

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_strains1);
hit_strains2 = translate(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

hit_strains1(inds) = [];
hit_data1(inds) = [];


inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

hit_strains2(inds) = [];
hit_data2(inds) = [];

hit_data1 = cell2mat(hit_data1);
hit_data2 = cell2mat(hit_data2);

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

[~,ind1,ind2] = intersect(hit_strains1, hit_strains2);
hit_strains_phaseI = hit_strains1(ind1);
hit_data_phaseI = hit_data1(ind1) - hit_data2(ind2);

%% Phase 2

[FILENAMES{end+1}, data3] = read_data('xlsread','./raw_data/Only original foldchanges/Direct comparison/Phase II/HOP_t0vsYPD10 (dir) PhII.xlsx', 'HOPt0vsYPD10_t0');
[FILENAMES{end+1}, data4] = read_data('xlsread','./raw_data/Only original foldchanges/Direct comparison/Phase II/HOPt0vsHOP10 (dir) PhII.xlsx', 'HOPt0vsHOP10_out0');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains3 = data3(2:end,1);
hit_data3 = data3(2:end,2);
essential3 = data3(2:end,5);

inds = find(strcmp('yes', essential3));
hit_strains3(inds) = [];
hit_data3(inds) = [];

hit_strains4 = data4(2:end,1);
hit_data4 = data4(2:end,3);
essential4 = data4(2:end,5);

inds = find(strcmp('yes', essential4));
hit_strains4(inds) = [];
hit_data4(inds) = [];
  
% Get ORFs
hit_strains3 = cellfun(@(x) regexp(x, '([\w\-])*(?=::)', 'match'), hit_strains3);
% hit_strains4 = cellfun(@(x) regexp(x, '([\w\-])*(?=::)', 'match'), hit_strains4);

% Eliminate all white spaces & capitalize
hit_strains3 = clean_orf(hit_strains3);
hit_strains4 = clean_orf(hit_strains4);

% If in gene name form, transform into ORF name
hit_strains3 = translate(hit_strains3);
hit_strains4 = translate(hit_strains4);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds));  

hit_strains3(inds) = [];
hit_data3(inds) = [];


inds = find(~is_orf(hit_strains4));
disp(hit_strains4(inds));  

hit_strains4(inds) = [];
hit_data4(inds) = [];

hit_data3 = cell2mat(hit_data3);
hit_data4 = cell2mat(hit_data4);

% If the same strain is present more than once, average its values
[hit_strains3, hit_data3] = grpstats(hit_data3, hit_strains3, {'gname','mean'});
[hit_strains4, hit_data4] = grpstats(hit_data4, hit_strains4, {'gname','mean'});

[~,ind1,ind2] = intersect(hit_strains3, hit_strains4);
hit_strains_phaseII = hit_strains3(ind1);
hit_data_phaseII = hit_data3(ind1) - hit_data4(ind2);

%%

hit_strains = unique([hit_strains_phaseI; hit_strains_phaseII]);
hit_data = nan(length(hit_strains), 2);

[~,ind1,ind2] = intersect(hit_strains_phaseI, hit_strains);
hit_data(ind2,1) = hit_data_phaseI(ind1);
[~,ind1,ind2] = intersect(hit_strains_phaseII, hit_strains);
hit_data(ind2,2) = hit_data_phaseII(ind1);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [197 198]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
novo_gonzalez_2013.orfs = hit_strains;
novo_gonzalez_2013.ph = hit_data_names;
novo_gonzalez_2013.data = hit_data;
novo_gonzalez_2013.dataset_ids = hit_data_ids;

%% Save

save('./novo_gonzalez_2013.mat','novo_gonzalez_2013');

%% Print out

fid = fopen('./novo_gonzalez_2013.txt','w');
write_matrix_file(fid, novo_gonzalez_2013.orfs, novo_gonzalez_2013.ph, novo_gonzalez_2013.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(novo_gonzalez_2013)
end

end

