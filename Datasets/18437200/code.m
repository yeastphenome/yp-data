%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
jin_freedman_2008.pmid = 18437200;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(jin_freedman_2008.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data1] = read_data('xlsread','./raw_data/journal.pgen.1000053.s005.xlsx', 'EC10');
[FILENAMES{end+1}, data2] = read_data('xlsread','./raw_data/journal.pgen.1000053.s005.xlsx', 'EC50');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_genenames1 = data1(5:end,2);
hit_genenames2 = data2(5:end,2);

% Get the data itself. Take the inverse and subtracting 1 because the original data
% represents GIF = (mut_unt/wt_unt)/(mut_treat/wt_treat) and higher values
% correspond to greater growth inhibition. Our convention is the opposite and expects no change to correspond to 0.

hit_data1 = 1./cell2mat(data1(5:end,3:9)) - 1;
hit_data2 = 1./cell2mat(data2(5:end,3:9)) - 1;

% Eliminate all white spaces & capitalize
hit_genenames1 = clean_genename(hit_genenames1);
hit_genenames2 = clean_genename(hit_genenames2);

hit_genenames2 = clean_orf(hit_genenames2); % to fix a few ORFs with "." instead of "-"

% If in gene name form, transform into ORF name
hit_strains1 = translate(hit_genenames1);
hit_strains2 = translate(hit_genenames2);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains1));
disp(hit_strains1(inds));  

inds = find(~is_orf(hit_strains2));
disp(hit_strains2(inds));  

% If the same strain is present more than once, average its values
[hit_strains1, hit_data1] = grpstats(hit_data1, hit_strains1, {'gname','mean'});
[hit_strains2, hit_data2] = grpstats(hit_data2, hit_strains2, {'gname','mean'});

% Merge 2 sets of data
hit_strains = unique([hit_strains1; hit_strains2]);
hit_data = zeros(length(hit_strains), size(hit_data1,2)+size(hit_data2,2));
[~,ind1,ind2] = intersect(hit_strains1, hit_strains);
hit_data(ind2,1:size(hit_data1,2)) = hit_data1(ind1,:);
[~,ind1,ind2] = intersect(hit_strains2, hit_strains);
hit_data(ind2,size(hit_data1,2)+1:end) = hit_data2(ind1,:);

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [11772 11773 11771 11774 11775 11776 11777 1311 1312 1313 1314 1310 1315 1316]';

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
jin_freedman_2008.orfs = hit_strains;
jin_freedman_2008.ph = hit_data_names;
jin_freedman_2008.data = hit_data;
jin_freedman_2008.dataset_ids = hit_data_ids;

%% Save

save('./jin_freedman_2008.mat','jin_freedman_2008');

%% Print out

fid = fopen('./jin_freedman_2008.txt','w');
write_matrix_file(fid, jin_freedman_2008.orfs, jin_freedman_2008.ph, jin_freedman_2008.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(jin_freedman_2008)
end

end

