%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
steinmetz_davis_2002.pmid = 12134146;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(steinmetz_davis_2002.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data (HOM)

[FILENAMES{end+1}, data1] = read_data('readtable','./raw_data/Regression_Tc1_hom.txt','Delimiter','\t','Format','%s %s %f %f %f %f %f %f %f %f %f %f');
[FILENAMES{end+1}, data2] = read_data('readtable','./raw_data/Regression_Tc2_hom.txt','Delimiter','\t','Format','%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f');

% Get the list of ORFs and the correponding data 
hit_strains1 = data1.orf;
hit_strains2 = data2.orf;

% Get the data itself
datasets1 = data1.Properties.VariableNames(8:end)';
datasets2 = data2.Properties.VariableNames(8:end)';
hit_data1 = data1{:,datasets1};
hit_data2 = data2{:,datasets2};

% Eliminate all white spaces & capitalize
hit_strains1 = clean_orf(hit_strains1);
hit_strains2 = clean_orf(hit_strains2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% Combine 1 and 2
hit_strains_hom = unique([hit_strains1; hit_strains2]);
datasets_hom = unique([datasets1; datasets2]);
hit_data_hom = nan(length(hit_strains_hom), length(datasets_hom), 2);

[~,ind1,ind2] = intersect(hit_strains1, hit_strains_hom);
[~,ind3,ind4] = intersect(datasets1, datasets_hom);
hit_data_hom(ind2,ind4,1) = hit_data1(ind1,ind3);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains_hom);
[~,ind3,ind4] = intersect(datasets2, datasets_hom);
hit_data_hom(ind2,ind4,2) = hit_data2(ind1,ind3);

hit_data_hom = nanmean(hit_data_hom,3);

% Remove rows that are all NaN
t = sum(~isnan(hit_data_hom),2);

hit_strains_hom(t==0) = [];
hit_data_hom(t==0,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids_hom = [4836 4840 4837 4838 4839 4841 4842 4844 4843]';

%% Load the data (HET)

[FILENAMES{end+1}, data3] = read_data('readtable','./raw_data/Regression_Tc1_het.txt','Delimiter','\t','Format','%s %s %f %f %f %f %f %f %f %f %f %f');
[FILENAMES{end+1}, data4] = read_data('readtable','./raw_data/Regression_Tc2_het.txt','Delimiter','\t','Format','%s %s %f %f %f %f %f %f %f %f %f %f %f');

% Get the list of ORFs and the correponding data 
hit_strains3 = data3.orf;
hit_strains4 = data4.orf;

% Get the data itself
datasets3 = data3.Properties.VariableNames(8:end)';
datasets4 = data4.Properties.VariableNames(8:end)';
hit_data3 = data3{:,datasets3};
hit_data4 = data4{:,datasets4};

% Eliminate all white spaces & capitalize
hit_strains3 = clean_orf(hit_strains3);
hit_strains4 = clean_orf(hit_strains4);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains3));
disp(hit_strains3(inds));  

inds = find(~is_orf(hit_strains4));
disp(hit_strains4(inds));  

% If the same strain is present more than once, average its values
[hit_strains3, hit_data3] = grpstats(hit_data3, hit_strains3, {'gname','mean'});
[hit_strains4, hit_data4] = grpstats(hit_data4, hit_strains4, {'gname','mean'});

% Combine 3 and 4
hit_strains_het = unique([hit_strains3; hit_strains4]);
datasets_het = unique([datasets3; datasets4]);
hit_data_het = nan(length(hit_strains_het), length(datasets_het), 2);

[~,ind1,ind2] = intersect(hit_strains3, hit_strains_het);
[~,ind3,ind4] = intersect(datasets3, datasets_het);
hit_data_het(ind2,ind4,1) = hit_data3(ind1,ind3);
[~,ind1,ind2] = intersect(hit_strains4, hit_strains_het);
[~,ind3,ind4] = intersect(datasets4, datasets_het);
hit_data_het(ind2,ind4,2) = hit_data4(ind1,ind3);

hit_data_het = nanmean(hit_data_het,3);

% Remove rows that are all NaN
t = sum(~isnan(hit_data_het),2);

hit_strains_het(t==0) = [];
hit_data_het(t==0,:) = [];

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids_het = [475 4831 4832 4830 4828 4827 4833 4834 4829]';

%% Prepare final dataset

hit_data_ids = [hit_data_ids_hom; hit_data_ids_het];

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

hit_strains = unique([hit_strains_hom; hit_strains_het]);
hit_data = nan(length(hit_strains), length(hit_data_ids));

[~,ind1,ind2] = intersect(hit_strains_hom, hit_strains);
hit_data(ind2,1:9) = hit_data_hom(ind1,:);
[~,ind1,ind2] = intersect(hit_strains_het, hit_strains);
hit_data(ind2,10:18) = hit_data_het(ind1,:);

% If the dataset is quantitative:
steinmetz_davis_2002.orfs = hit_strains;
steinmetz_davis_2002.ph = hit_data_names;
steinmetz_davis_2002.data = hit_data;
steinmetz_davis_2002.dataset_ids = hit_data_ids;

%% Save

save('./steinmetz_davis_2002.mat','steinmetz_davis_2002');

%% Print out

fid = fopen('./steinmetz_davis_2002.txt','w');
write_matrix_file(fid, steinmetz_davis_2002.orfs, steinmetz_davis_2002.ph, steinmetz_davis_2002.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(steinmetz_davis_2002)
end

end

