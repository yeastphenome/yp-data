%% McCormick~Kennedy, 2015

function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
mccormick_kennedy_2015.pmid = 26456335;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(mccormick_kennedy_2015.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

[FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/rls-summary-for-Anastasia-Baryshnikova-all-BY-haploid-deletion-YPD-30C-mm042018.xlsx', 'rls');

% Eliminate all strains that are not BY4742 or BY4741
data(cellfun(@isnumeric, data(:,5)),5) = {''};
data(strcmp(data(:,5), 'BY4,742'),5) = {'BY4742'};
data = data(ismember(data(:,5),{'BY4741','BY4742'}),:);

data(cellfun(@isnumeric, data(:,6)),6) = {''};
data(strcmp(data(:,6), 'matalpha'),6) = {'MATalpha'};
data = data(ismember(data(:,6),{'MATa','MATalpha'}),:);

% Eliminate everything BUT the single mutants
t = cellfun(@strtrim, data(:,8), 'UniformOutput', false);
t = cellfun(@(x) ~contains(x, ' '), t, 'UniformOutput', false);
t = cell2mat(t);
data = data(t,:);

% Separate the 2 backgrounds
data_a = data(strcmp(data(:,5), 'BY4741'),:);
data_b = data(strcmp(data(:,5), 'BY4742'),:);

% Only keep the BY4742 background (systematic screen)

% Merge the replicates
data_b_gene = clean_genename(data_b(:,8));
data_b_orf = translate(data_b_gene);

inds = find(~is_orf(data_b_orf));
data_b_orf(inds) = [];
data_b(inds,:) = [];

data_b_u_orf = unique(data_b_orf);
data_b_u_genename = translate(data_b_u_orf, 'genenames');
data_b_u_rls = cell(length(data_b_u_orf),1);
data_b_u_rls_ref = cell(length(data_b_u_orf),1);

for i = 1 : length(data_b_u_orf)
    data_b_u_rls{i} = [];
    inds = find(strcmp(data_b_u_orf{i}, data_b_orf));
    for j = 1 : length(inds)
        if ~isnumeric(data_b{inds(j),15})
            lifespans = cellfun(@str2num, strsplit(data_b{inds(j), 15}, ','));
        else
            lifespans = [data_b{inds(j), 15}];
        end
        data_b_u_rls{i} = cat(2, data_b_u_rls{i}, lifespans);
        
        % Ref lifespans
        if ~isnumeric(data_b{inds(j),28})
            lifespans = cellfun(@str2num, strsplit(data_b{inds(j), 28}, ','));
        else
            lifespans = [data_b{inds(j), 28}];
        end
        data_b_u_rls_ref{i} = cat(2, data_b_u_rls_ref{i}, lifespans);
        
    end
end

data_b_u_rls_mean = cellfun(@nanmean, data_b_u_rls);
data_b_u_rls_std = cellfun(@nanstd, data_b_u_rls);
data_b_u_rls_nr = cellfun(@numel, data_b_u_rls);

data_b_u_rls_ref_mean = cellfun(@nanmean, data_b_u_rls_ref);
data_b_u_rls_ref_std = cellfun(@nanstd, data_b_u_rls_ref);
data_b_u_rls_ref_nr = cellfun(@numel, data_b_u_rls_ref);

% bins = [15:1:40];
% h1 = histc(data_b_u_rls_mean, bins);
% h2 = histc(data_b_u_rls_ref_mean, bins);
% 
% figure()
% bar(bins, h1/sum(h1),'FaceColor','b','FaceAlpha', 0.5);
% hold all;
% bar(bins, h2/sum(h2),'FaceColor','r', 'FaceAlpha', 0.5);

% Normalize to WT
data_b_u_rls_ratio = data_b_u_rls_mean ./ data_b_u_rls_ref_mean;
data_b_u_rls_ratio_std = data_b_u_rls_ratio .* sqrt((data_b_u_rls_std./data_b_u_rls_mean).^2 + (data_b_u_rls_ref_std./data_b_u_rls_ref_mean).^2);

% Only keep data with n > 5 (rest is unreliable)
% The authors subsequently retested these strains but the raw version of that data is
% (unfortunately) not recoverable. Only the published results.
% Since all but one of the published (most reliable) strains has n > 5,
% we've decided to set the n < 5 strains to 0 (instead of NaN) to indicate that they are
% likely neither short-lived nor long-lived (instead of "not tested")

inds = find(data_b_u_rls_nr <= 5 & data_b_u_rls_ref_nr <= 5);
data_b_u_rls_ratio(inds) = 0;
data_b_u_rls_ratio_std(inds) = 0;

data_b_u_rls_mean(inds) = 0;
data_b_u_rls_std(inds) = 0;
data_b_u_rls_nr(inds) = 0;

data_b_u_rls_ref_mean(inds) = 0;
data_b_u_rls_ref_std(inds) = 0;
data_b_u_rls_ref_nr(inds) = 0;

%% The published data is effectively a subset of the raw data
% 
% [FILENAMES{end+1}, data_hits] = read_data('xlsread','./raw_data/mmc3.xlsx', 'Table S2');
% 
% % Get the list of ORFs and the correponding data 
% hit_strains = data_hits(4:end,1);
% 
% % Get the data itself
% hit_data = cell2mat(data_hits(4:end, 12));
%    
% % Eliminate all white spaces & capitalize
% hit_strains = clean_genename(hit_strains);
% 
% % If in gene name form, transform into ORF name
% hit_strains = translate(hit_strains);
% 
% % Find anything that doesn't look like an ORF
% hit_strains(ismember(hit_strains, {'FMP42'})) = {'YMR221C'};
% inds = find(~is_orf(hit_strains));
% hit_strains(inds) = [];
% hit_data(inds) = [];

%% Prepare dataset

hit_strains = data_b_u_orf;
hit_data = data_b_u_rls_ratio;

% MANUAL. Get the dataset ids corresponding to each dataset (in order)
% Multiple datasets (e.g., replicates) may get the same id, which can then
% be used to average them out
hit_data_ids = [696];

%% Tested strains (as published)
% 
% % Load tested strains
% [FILENAMES{end+1}, data] = read_data('xlsread','./raw_data/mmc2.xlsx', 'Table S1');
% 
% % Get the list of ORFs and the correponding data 
% tested_strains = data(3:end,1);
% 
% % Eliminate all white spaces & capitalize
% tested_strains = clean_orf(tested_strains);
% 
% % Find anything that doesn't look like an ORF
% inds = find(~is_orf(tested_strains));
% tested_strains(inds) = [];
% 
% % Finally, take the unique set
% tested_strains = unique(tested_strains);
% 
% % Make sure the that all the hits are part of the tested set
% [missing,~] = setdiff(hit_strains, tested_strains);
% 
% % If it seems reasonable, add the missing hits to the list of tested strains
% tested_strains = [tested_strains; missing];

%% Prepare final dataset

% Match the dataset ids with the dataset standard names
[~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
hit_data_names = cell(size(hit_data_ids));
hit_data_names(ind2) = datasets.standard_name(ind1);

% If the dataset is discrete/binary and the tested strains were provided separately:
mccormick_kennedy_2015.orfs = hit_strains;
mccormick_kennedy_2015.ph = hit_data_names;
mccormick_kennedy_2015.data = hit_data;
mccormick_kennedy_2015.dataset_ids = hit_data_ids;

%% Save

save('./mccormick_kennedy_2015.mat','mccormick_kennedy_2015');

%% Print out

fid = fopen('./mccormick_kennedy_2015.txt','w');
write_matrix_file(fid, mccormick_kennedy_2015.orfs, mccormick_kennedy_2015.ph, mccormick_kennedy_2015.data);
fclose(fid);

%% Save to DB (admin)

addpath(genpath('../../Private-Utils/'));
if exist('save_data_to_db.m')
    res = save_data_to_db(mccormick_kennedy_2015)
end

end

