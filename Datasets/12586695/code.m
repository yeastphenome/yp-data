%% FirstAuthor~LastAuthor, YYYY
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
enyenihi_saunders_2003.pmid = 12586695;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(enyenihi_saunders_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/Table_1_of_total_results.xlsx', 'Figure 1');

% Get the list of ORFs and the correponding data 
% (this part usually changes significantly based on the format of the raw data file)
hit_strains = data(:,1:2:31);
hit_strains = hit_strains(:);

% Get the data itself
hit_data = [data(:,2:2:32)];
hit_data = hit_data(:);
   
% Eliminate all white spaces & capitalize
inds = find(cellfun(@isnumeric, hit_strains));
hit_strains(inds) = [];
hit_data(inds) = [];

hit_strains = clean_genename(hit_strains);

% If in gene name form, transform into ORF name
hit_strains = translate(hit_strains);

% If possible, fix the problem (typos, omissions etc.)
hit_strains(ismember(hit_strains, {'YBLO29C'})) = {'YBL029C'};
hit_strains(ismember(hit_strains, {'KRE20'})) = {'YAL056C-A'};
hit_strains(ismember(hit_strains, {'KRE23'})) = {'YAL042C-A'};
hit_strains(ismember(hit_strains, {'KRE24'})) = {'YPL102C'};
hit_strains(ismember(hit_strains, {'FYV14'})) = {'YDL213C'};
hit_strains(ismember(hit_strains, {'FYV9'})) = {'YDR140W'};
hit_strains(ismember(hit_strains, {'FYV11'})) = {'YFL023W'};
hit_strains(ismember(hit_strains, {'FYV2'})) = {'YIL054W'};
hit_strains(ismember(hit_strains, {'MFALPHA2'})) = {'YGL089C'};
hit_strains(ismember(hit_strains, {'RPSOA'})) = {'YGR214W'};
hit_strains(ismember(hit_strains, {'HSP-150'})) = {'YJL159W'};
hit_strains(ismember(hit_strains, {'PET1OO'})) = {'YDR079W'};
hit_strains(ismember(hit_strains, {'EFR4'})) = {'YLR114C'};
hit_strains(ismember(hit_strains, {'WH12'})) = {'YOR043W'};
hit_strains(ismember(hit_strains, {'MFALPHA1'})) = {'YPL187W'};
hit_strains(ismember(hit_strains, {'SHD7'})) = {'YPL180W'};
hit_strains(ismember(hit_strains, {'SDF1'})) = {'YPR040W'};
hit_strains(find(strncmp('YDL062W', hit_strains, 7))) = {'YDL062W'};

% Find anything that doesn't look like an ORF
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

hit_strains(inds) = [];
hit_data(inds) = [];

hit_data(strcmp('L0W', hit_data)) = {'LOW'};
hit_data = clean_genename(hit_data);    % a quick way to remove unwanted spaces etc.

hit_data(strcmp('normal', hit_data)) = {0};
hit_data(strncmp('LOW', hit_data, 3)) = {-1};
hit_data(strncmp('VLOW', hit_data, 4)) = {-2};
hit_data(strcmp('NONE', hit_data)) = {-3};

inds = find(~cellfun(@isnumeric, hit_data));
hit_strains(inds) = [];
hit_data(inds) = [];

hit_data = cell2mat(hit_data);

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [72];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
enyenihi_saunders_2003.orfs = hit_strains;
enyenihi_saunders_2003.ph = hit_data_names;
enyenihi_saunders_2003.data = hit_data;
enyenihi_saunders_2003.dataset_ids = hit_data_ids;

%% Save

save('./enyenihi_saunders_2003.mat','enyenihi_saunders_2003');

%% Print out

fid = fopen('./enyenihi_saunders_2003.txt','w');
write_matrix_file(fid, enyenihi_saunders_2003.orfs, enyenihi_saunders_2003.ph, enyenihi_saunders_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(enyenihi_saunders_2003)
end

end

