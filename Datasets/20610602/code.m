%% Cooper~Fields, 2010
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

cooper_fields_2010.pmid = 20610602;

[FILENAMES{end+1}, data.raw] = dataread('xlsread','./raw_data/SupplementalTable4.xlsx');

cooper_fields_2010.orfs = upper(data.raw(2:end,1));
cooper_fields_2010.ph = data.raw(1,3:end)';

inds = find(strcmp('-', data.raw));
data.raw(inds) = {NaN};

cooper_fields_2010.data = cell2mat(data.raw(2:end,3:end));
cooper_fields_2010.desc = {'Supplemental Table 4. Log2 transformed ratios representing fold change compared to average for each identified amino acid in each of nearly 4500 samples.'};


% Identical ORFs with multiple values -> average
[t,t2] = grpstats(cooper_fields_2010.data, cooper_fields_2010.orfs, {'gname','mean'});
cooper_fields_2010.orfs = t;
cooper_fields_2010.data = t2;

a = mfilename('fullpath');
a = a(1:end-4);
save([a,'cooper_fields_2010.mat'],'cooper_fields_2010');
return;

% Save data into database
insert_data_into_db(cooper_fields_2010);

end

