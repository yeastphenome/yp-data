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
hap_hits = {};
hom_hits = {};
all_het_hits = {};
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            hit_strains = [hit_strains; C{j}];
            each_data = regexp(data{i}, '(wt|\d+)\s(wt|\d+)\s(wt|\d+|NA)', 'tokens');
            each_data = each_data{1};
            for k = 1:3
                if ~isempty(strfind(each_data{k}, 'wt'))
                    each_data(k) = {'100'};
                end
                if ~isempty(strfind(each_data{k}, 'NA'))
                    each_data(k) = {'NaN'};
                end
                if (k == 1)
                    hap_hits = [hap_hits; cell2mat(each_data(k))];
                end
                if (k == 2)
                    hom_hits = [hom_hits; cell2mat(each_data(k))];
                end
                if (k == 3)
                    all_het_hits = [all_het_hits; cell2mat(each_data(k))];
                end
            end
        end
    end
end


%% Load tested - Haploinsufficiency

[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/Table 4.txt', '%s', 'Delimiter', '\n');

% Resistant
data_r = data(5:37);
het_hits = {};
het_strains = {};
for i = 1:length(data_r)
    C = strsplit(data_r{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            het_strains = [het_strains; C{j}];
            each_data = regexp(data_r{i}, '\s(\d\d+)', 'tokens');
            each_data = each_data{1};
            het_hits = [het_hits; cell2mat(each_data(1))];
        end
    end
end

% Sensitive
data_s = data(38:end);
for i = 1:length(data_s)
    C = strsplit(data_s{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            het_strains = [het_strains; C{j}];
            each_data = regexp(data_r{i}, '\s(\d\d+)', 'tokens');
            each_data = each_data{1};
            het_hits = [het_hits; cell2mat(each_data(1))];
        end
    end
end

% All in het_strains are unique.

%% Load tested - Hypersensitive

[FILENAMES{end+1}, data] = read_data('textscan', './raw_data/Table 3.txt', '%s', 'Delimiter', '\n');

%Isolate the ORFs from rest of the data
for i = 1:length(data)
    C = strsplit(data{i}, ' ');
    for j = 1:length(C)
        if is_orf(C(j))
            hit_strains = [hit_strains; C{j}];
            each_data = regexp(data{i}, '(wt|\d+)\s(wt|\d+)\s(wt|\d+|NA)', 'tokens');
            each_data = each_data{1};
            for k = 1:3
                if ~isempty(strfind(each_data{k}, 'wt'))
                    each_data(k) = {'100'};;
                end
                if ~isempty(strfind(each_data{k}, 'NA'))
                    each_data(k) = {'NaN'};;
                end
                if (k == 1)
                    hap_hits = [hap_hits; cell2mat(each_data(k))];
                end
                if (k == 2)
                    hom_hits = [hom_hits; cell2mat(each_data(k))];
                end
                if (k == 3)
                    all_het_hits = [all_het_hits; cell2mat(each_data(k))];
                end
            end
        end
    end
end

%% Add the data together

% Convert cell arrays to numeric vectors
all_het_hits = cell2mat(all_het_hits);
hap_hits = cell2mat(hap_hits);
hom_hits = cell2mat(hom_hits);

% add het_strains to hit_strains
% then add het_hits to all_het_hits
hit_strains = [hit_strains; het_strains];

%% Here is where I need to put the data in, but I'm not sure how to 
%% do that with our several collections


%% Prepare Final Dataset

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [81, 82, 83];

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

end
