%% Cooper~Fields, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

cooper_fields_2010.pmid = 20610602;

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/SupplementalTable4.xlsx');

cooper_fields_2010.orfs = data.raw(2:end,1);
cooper_fields_2010.orfs = clean_orf(cooper_fields_2010.orfs);

inds = find(~is_orf(cooper_fields_2010.orfs));
disp(cooper_fields_2010.orfs(inds));

% Fix the typo
cooper_fields_2010.orfs(strcmp('YML048WA-', cooper_fields_2010.orfs)) = {'YML048W-A'};

cooper_fields_2010.ph = data.raw(1,3:end)';

inds = find(strcmp('-', data.raw));
data.raw(inds) = {NaN};

cooper_fields_2010.data = cell2mat(data.raw(2:end,3:end));
cooper_fields_2010.desc = {'Supplemental Table 4. Log2 transformed ratios representing fold change compared to average for each identified amino acid in each of nearly 4500 samples.'};

% Load the phenotype mapping file
[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/phenotype_mapping.xlsx','Sheet1');

ph.orig = data(:,1);
ph.who = data(:,2);
ph.what = data(:,3);
ph.where = data(:,4);
ph.when = data(:,5);
ph.how = data(:,6);

ph = build_phenotype_name(ph);

[~,ind1,ind2] = intersect(cooper_fields_2010.ph, ph.orig);
cooper_fields_2010.ph(ind1) = ph.ph(ind2);

% Average the data for the 2 lysine peaks
[t, t2] = grpstats(cooper_fields_2010.data', cooper_fields_2010.ph, {'gname','mean'});
cooper_fields_2010.ph = t;
cooper_fields_2010.data = t2';

% Identical ORFs with multiple values -> average
[t,t2] = grpstats(cooper_fields_2010.data, cooper_fields_2010.orfs, {'gname','mean'});
cooper_fields_2010.orfs = t;
cooper_fields_2010.data = t2;

save('./cooper_fields_2010.mat','cooper_fields_2010');

end

