%% Chan~Zheng, 2000
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};

chan_zheng_2000.source = {'Supplementary Table 3'};
chan_zheng_2000.downloaddate = {'2014-01-31'};
chan_zheng_2000.pmid = 11078525;

[FILENAMES{end+1}, data.raw] = read_data('xlsread','./raw_data/chan_zheng_2000_HAP.xlsx', 'Sheet1');

phenotypes = {'growth'};
treatments = {'rapamycin, 25 nM'};

% Eliminate all white spaces & capitalize
data.raw(:,1) = upper(clean_orf(data.raw(:,1)));

chan_zheng_2000.orfs = data.raw(:,1);
chan_zheng_2000.data = cell2mat(data.raw(:,2));
chan_zheng_2000.ph = [strcat(phenotypes, '; ', treatments)];

save('./chan_zheng_2000.mat','chan_zheng_2000');
return;

% Save data into database
dt = chan_zheng_2000;
datasets = get_datasets_for_paper(dt);

[~,database_ix] = sortrows(datasets.names,[1 2 3]);
[~,ph_ix] = sort(dt.ph);

% Before loading into database, manually check the order of ph_ix and database_ix to make sure they correspond.
datasets.names(database_ix,:)
dt.ph(ph_ix)

insert_data_into_db(dt, ph_ix, datasets.ids(database_ix));

fid = fopen('./chan_zheng_2000.txt','w');
write_matrix_file(fid, chan_zheng_2000.orfs, chan_zheng_2000.ph, chan_zheng_2000.data);
fclose(fid);

end

