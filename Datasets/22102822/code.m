%% Berry~Gasch, 2011
function FILENAMES = code()
addpath(genpath('../../Yeast-Matlab-Utils/'));
FILENAMES = {};
berry_gasch_2011.pmid = 22102822;

phenotypes = {'growth'};
treatment = {'control', '0.4 mM H202', '2.5 mM DDT', '0.7 M NaCl', 'Heat Shock', '20 uM tunicamycin', '2.5 mM DDT, H2O2', '0.7 M NaCl, H2O2', 'Heat Shock, H2O2', '20 uM tunicamycin, H2O2'};

%% Hit Strains

% Load file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/pgen.1002353.s009.xlsx', 'Hom-Het COMPILATION');

% Get the list of ORFs
strains = data(2:end, 1);

% Clean up ORFs
strains = clean_orf(strains);

% Find anything that doesn't look like an ORF
inds = find(~is_orf(strains));
disp(strains(inds)); 

% Get data from hits
hit_data = data(2:end, 2:end);
hit_data(strcmp('NA', hit_data)) = {NaN};
indx = ~cellfun(@isnumeric, hit_data);
hit_data(indx) = {NaN};
hit_data = cell2mat(hit_data);

%% Get Conditions and Corresponding data
hit_cond = data(1, 2:end);

% Condition 1: Sample 1 vs Sample 0
% find all the conditions with "sample0"
indx = strfind(hit_cond, 'Sample0');
indx = find(~(cellfun(@isempty, indx)));
cond_one_data = nanmean(hit_data(:,indx),2);

% Condition 2: Sample 2 vs Sample 1
% find all the conditions with "sample2"
indx = strfind(hit_cond, 'Sample2');
indx = find(~(cellfun(@isempty, indx)));
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

% D. the conditions with TM
indx = find(~cellfun(@isempty, regexp(hit_cond, 'TM[0-9] Sample4[A-Z]? vs Sample3 Array')));
cond_fourD_data = nanmean(hit_data(:,indx),2);

final_data = [cond_one_data cond_two_data cond_threeA_data cond_threeB_data cond_threeC_data cond_threeD_data cond_fourA_data cond_fourB_data cond_fourC_data cond_fourD_data];

% Average any repeated value
[strains, final_data] = grpstats(final_data, strains, {'gname','mean'});

% Prepare final dataset
berry_gasch_2011.orfs = strains;
berry_gasch_2011.ph = strcat(phenotypes, '; ', treatment);
berry_gasch_2011.data = final_data;

%% Save

save('./berry_gasch_2011.mat','berry_gasch_2011');

fid = fopen('./berry_gasch_2011.txt','w');
write_matrix_file(fid, berry_gasch_2011.orfs, berry_gasch_2011.ph, berry_gasch_2011.data);
fclose(fid);

end