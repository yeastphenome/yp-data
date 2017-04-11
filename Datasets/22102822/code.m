%% Berry~Gasch, 2011
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
berry_gasch_2011.pmid = 22102822;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(berry_gasch_2011.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pgen.1002353.s009.xlsx', 'Hom-Het COMPILATION');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = clean_orf(strains);

strains = translate(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = data(2:end, 2:end);
hit_data(strcmp('NA', hit_data)) = {NaN};

indx = find(~cellfun(@isnumeric, hit_data));
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

%% Get Conditions and Corresponding data

hit_cond = data(1, 2:end);

% Sample0 = unstressed, time 0 control sample
% Sample1 = outgrown in YPD for 10 generation (with dilutions to maintain log phase growth), compared to Sample0 to identify slow-growing strains
% Sample2 = exposed to 0.4 mM H2O2 for 2 h, outgrown for 10 generations
% Sample3 = exposed to 1 of 3 primary stresses (NaCl, DTT, 40C), outgrown for 10 generations
% Sample4 = exposed to secondary stress (1 mM or 1.2 mM H2O2), outgrown for 10 generations

% MicroArray Data
% Condition 1: Sample 1 vs Sample 0
% find all the conditions with "sample0"
indx = find(~cellfun(@isempty, regexp(hit_cond, 'Sample1 vs Sample0 Array')));
cond_one_data = nanmean(hit_data(:,indx),2);

% Condition 2: Sample 2 vs Sample 1
% find all the conditions with "sample2"
indx = find(~cellfun(@isempty, regexp(hit_cond, 'Sample2 vs Sample1 Array')));
cond_two_data = nanmean(hit_data(:,indx),2);

% Condition 3: Sample 3 vs Sample 1
% should have 4 different sets of these
% A. the conditions with DTT
indx = find(~cellfun(@isempty, regexp(hit_cond, 'DTT[0-9] Sample3 vs Sample1 Array')));
cond_threeA_data = nanmean(hit_data(:,indx),2);

% B. the conditions with NaCl
indx = find(~cellfun(@isempty, regexp(hit_cond, 'NaCl[0-9] Sample3 vs Sample1 Array')));
cond_threeB_data = nanmean(hit_data(:,indx),2);

% C. the conditions with HS
indx = find(~cellfun(@isempty, regexp(hit_cond, 'HS[0-9] Sample3 vs Sample1 Array')));
cond_threeC_data = nanmean(hit_data(:,indx),2);

% D. the conditions with TM
indx = find(~cellfun(@isempty, regexp(hit_cond, 'TM[0-9] Sample3 vs Sample1 Array')));
cond_threeD_data = nanmean(hit_data(:,indx),2);

% Condition 4: Sample 4 vs Sample 3
% should have 4 different sets of these
% A. the conditions with DTT
indx = find(~cellfun(@isempty, regexp(hit_cond, 'DTT[0-9] Sample4[A-Z]? vs Sample3 Array')));
cond_fourA_data = nanmean(hit_data(:,indx),2);

% B. the conditions with NaCl
indx = find(~cellfun(@isempty, regexp(hit_cond, 'NaCl[0-9] Sample4[A-Z]? vs Sample3 Array')));
cond_fourB_data = nanmean(hit_data(:,indx),2);

% C. the conditions with HS
indx = find(~cellfun(@isempty, regexp(hit_cond, 'HS[0-9] Sample4[A-Z]? vs Sample3 Array')));
cond_fourC_data = nanmean(hit_data(:,indx),2);

% % D. the conditions with TM
% indx = find(~cellfun(@isempty, regexp(hit_cond, 'TM[0-9] Sample4[A-Z]? vs Sample3 Array')));        % <<<< EMPTY
% cond_fourD_data = nanmean(hit_data(:,indx),2);

% Sequencing
% Condition 1: Sample 1 vs Sample 0
% find all the conditions with "sample0"
indx = find(~cellfun(@isempty, regexp(hit_cond, 'Sample1 vs Sample0 [(DN)(UP)]')));
cond_one_data_seq = nanmean(hit_data(:,indx),2);

% Condition 2: Sample 2 vs Sample 1
% find all the conditions with "sample2"
indx = find(~cellfun(@isempty, regexp(hit_cond, 'Sample2 vs Sample1 [(DN)(UP)]')));
cond_two_data_seq = nanmean(hit_data(:,indx),2);

% Condition 3: Sample 3 vs Sample 1
% should have 4 different sets of these
% A. the conditions with DTT
indx = find(~cellfun(@isempty, regexp(hit_cond, 'DTT[0-9] Sample3 vs Sample1 [(DN)(UP)]')));
cond_threeA_data_seq = nanmean(hit_data(:,indx),2);

% B. the conditions with NaCl
indx = find(~cellfun(@isempty, regexp(hit_cond, 'NaCl[0-9] Sample3 vs Sample1 [(DN)(UP)]')));
cond_threeB_data_seq = nanmean(hit_data(:,indx),2);

% C. the conditions with HS
indx = find(~cellfun(@isempty, regexp(hit_cond, 'HS[0-9] Sample3 vs Sample1 [(DN)(UP)]')));
cond_threeC_data_seq = nanmean(hit_data(:,indx),2);

% D. the conditions with TM
indx = find(~cellfun(@isempty, regexp(hit_cond, 'TM[0-9] Sample3 vs Sample1 [(DN)(UP)]')));
cond_threeD_data_seq = nanmean(hit_data(:,indx),2);

% Condition 4: Sample 4 vs Sample 3
% should have 4 different sets of these
% A. the conditions with DTT
indx = find(~cellfun(@isempty, regexp(hit_cond, 'DTT[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]')));
cond_fourA_data_seq = nanmean(hit_data(:,indx),2);

% B. the conditions with NaCl
indx = find(~cellfun(@isempty, regexp(hit_cond, 'NaCl[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]')));
cond_fourB_data_seq = nanmean(hit_data(:,indx),2);

% C. the conditions with HS
indx = find(~cellfun(@isempty, regexp(hit_cond, 'HS[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]')));
cond_fourC_data_seq = nanmean(hit_data(:,indx),2);

% D. the conditions with TM
indx = find(~cellfun(@isempty, regexp(hit_cond, 'TM[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]')));
cond_fourD_data_seq = nanmean(hit_data(:,indx),2);


final_data = [cond_one_data cond_two_data ...
              cond_threeA_data cond_threeB_data cond_threeC_data cond_threeD_data ...
              cond_fourA_data cond_fourB_data cond_fourC_data ...
              cond_one_data_seq cond_two_data_seq ...
              cond_threeA_data_seq cond_threeB_data_seq cond_threeC_data_seq cond_threeD_data_seq ...
              cond_fourA_data_seq cond_fourB_data_seq ...
              cond_fourC_data_seq cond_fourD_data_seq];
          
hit_data_ids = [758; 759; 761; 760; 762; 763; 590; 588; 589; 5395; 5396; 5397; 5398; 5399; 5400; 5401; 5402; 5403; 5404];

% Average any repeated value
[strains, final_data] = grpstats(final_data, strains, {'gname','mean'});

% To separate HOM from HET data, identify the essential genes and remove
% them (for now)

load essential_genes_151215.mat
[~,ind1,ind2] = intersect(essential_genes, strains);
strains(ind2) = [];
final_data(ind2,:) = [];


%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is quantitative:
berry_gasch_2011.orfs = strains;
berry_gasch_2011.ph = hit_data_names;
berry_gasch_2011.data = final_data;
berry_gasch_2011.dataset_ids = hit_data_ids;

%% Save

save('./berry_gasch_2011.mat','berry_gasch_2011');

%% Print out

fid = fopen('./berry_gasch_2011.txt','w');
write_matrix_file(fid, berry_gasch_2011.orfs, berry_gasch_2011.ph, berry_gasch_2011.data);
fclose(fid);

end
