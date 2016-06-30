%% Page~Bussey, 2003
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
page_bussey_2003.pmid = 12663529;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(page_bussey_2003.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load tested - Resistant

[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/Table 2.txt', '%s', 'Delimiter', '\n');

%Isolate the ORFs from rest of the data
hit_strains = {};
hap_data = [];
hom_data = [];
het_data = [];
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            hit_strains = [hit_strains; C{j}];
            each_data = regexp(data{i}, '(wt|\d+)\s(wt|\d+)\s(wt|\d+|NA)', 'tokens');
            each_data = each_data{1};
            for k = 1:3
                if strcmp(each_data{k}, 'wt')
                    each_data(k) = {'100'};
                end
                if strcmp(each_data{k}, 'NA')
                    each_data(k) = {'NaN'};
                end
                if k == 1
                    hap_data = [hap_data; str2num(each_data{k})];
                end
                if k == 3
                    hom_data = [hom_data; str2num(each_data{k})];
                end
                if k == 2
                    het_data = [het_data; str2num(each_data{k})];
                end
            end
        end
    end
end


%% Load tested - Haploinsufficiency

[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/Table 4.txt', '%s', 'Delimiter', '\n');

% Resistant
data = data(5:end);
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            hit_strains = [hit_strains; C{j}];
            each_data = regexp(data{i}, '\s(wt|\d+)', 'tokens');
            each_data = each_data{1};
            if strcmp(each_data{1}, 'wt')
                each_data(1) = {'100'};
            end
            het_data = [het_data; str2num(each_data{1})];
            hap_data = [hap_data; NaN];
            hom_data = [hom_data; NaN];
        end
    end
end

%% Load tested - Hypersensitive

[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/Table 3.txt', '%s', 'Delimiter', '\n');

%Isolate the ORFs from rest of the data
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            hit_strains = [hit_strains; C{j}];
            each_data = regexp(data{i}, '(wt|\s\d+)\s(wt|\d+)\s(wt|\d+|NA)', 'tokens');
            each_data = each_data{1};
            for k = 1:3
                if strcmp(each_data{k}, 'wt')
                    each_data(k) = {'100'};
                end
                if strcmp(each_data{k}, 'NA')
                    each_data(k) = {'NaN'};
                end
                if k == 1
                    hap_data = [hap_data; str2num(each_data{k})];
                end
                if k == 3
                    hom_data = [hom_data; str2num(each_data{k})];
                end
                if k == 2
                    het_data = [het_data; str2num(each_data{k})];
                end
            end
        end
    end
end


%% Do the regular checks

hit_strains = clean_orf(hit_strains);
inds = find(~is_orf(hit_strains));
disp(hit_strains(inds));

hit_data = [hap_data het_data hom_data];

% If the same strain is present more than once, average its values
[hit_strains, hit_data] = grpstats(hit_data, hit_strains, {'gname','mean'});

% Transform the data, so that resistant strains have positive values, 
% sensitive strains have negative values, and wt = 0
hit_data = 100 - hit_data;

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [81; 83; 82];

%% Prepare Final Dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
page_bussey_2003.orfs = hit_strains;
page_bussey_2003.ph = hit_data_names;
page_bussey_2003.data = hit_data;
page_bussey_2003.dataset_ids = hit_data_ids;

%% Save

save('./page_bussey_2003.mat','page_bussey_2003');

%% Print out

fid = fopen('./page_bussey_2003.txt','w');
write_matrix_file(fid, page_bussey_2003.orfs, page_bussey_2003.ph, page_bussey_2003.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(page_bussey_2003)
end

end
