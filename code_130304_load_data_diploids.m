% %% Enyenihi~Saunders, 2003
% % DATA = enyenihi_saunders_2003
% 
% enyenihi_saunders_2003.source = {'http://www.genetics.org/content/suppl/2003/02/21/163.1.47.DC1/Table_1_of_total_results.xls'};
% enyenihi_saunders_2003.downloaddate = {'2014-02-13'};
% enyenihi_saunders_2003.pmid = 12586695;
% 
% [data.num,data.txt,data.raw] = ...
%     xlsread('Datasets/Phenotypes/2003_Enyenihi~Saunders/Table_1_of_total_results.xlsx','Figure 1');
% 
% phenotypes = {'sporulation'};
% treatments = {''};
% 
% [old_genename, new_orf] = textread('Datasets/Phenotypes/2003_Enyenihi~Saunders/genename_updates_SGD140214.txt','%s %s');
% for chr = 1 : 16
%     genenames = data.raw(2:end,chr*2-1);
%     results = data.raw(2:end,chr*2);
%     inds = find(cellfun(@isnumeric, genenames));
%     genenames(inds) = [];
%     results(inds) = [];
%     genenames = cellfun(@strtrim, genenames,'UniformOutput',0);
%     inds = find(strcmp('centromere',lower(genenames)));
%     genenames(inds) = [];
%     results(inds) = [];
%     [t,ind1,ind2] = intersect(genenames, old_genename);
%     genenames(ind1) = new_orf(ind2);
%     
%     [orfs, unmatched] = genename2orf(genenames);
%     
% end
% 
% % Eliminate white spaces before/after ORF
% data.raw(2:end,1) = cellfun(@strtrim, data.raw(2:end,1),'UniformOutput',0);
% 
% % Replace 'Inf' with Inf
% t = data.raw(2:end,2:end);
% t(~cellfun(@isnumeric, t)) = {Inf};
% 
% blackburn_avery_2003.orfs = upper(data.raw(2:end,1));
% blackburn_avery_2003.data = cell2mat(t);
% 
% blackburn_avery_2003.ph = [strcat(phenotypes{1}, '; ', treatments)];
% 
% 
% % Eliminate the essential genes
% load essential_genes_100908;
% [t,ind1,ind2] = intersect(blackburn_avery_2003.orfs, essential_genes);
% blackburn_avery_2003.orfs(ind1) = [];
% blackburn_avery_2003.data(ind1,:) = [];
% 
% save('Datasets/Phenotypes/2003_Blackburn~Avery/blackburn_avery_2003.mat','blackburn_avery_2003'); 


%% Hoepfner~Movva, 2014
% DATA = hoepfner_movva_2014
hoepfner_movva_2014.pmid = 24360837;

hop_data = read_matrix_file('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Hoepfner~Movva/HOP_scores.txt',1,1);

hop_data.labels_col = cellfun(@(x) regexprep(x,'"',''), hop_data.labels_col, 'UniformOutput',0);    % remove the unnecessary quotes from strings

% Extract coumpound ID and concentration
hop_data.cmd = nan(size(hop_data.labels_col));
hop_data.conc = nan(size(hop_data.labels_col));
hop_data.group = cell(size(hop_data.labels_col));
hop_data.zscore = cell(size(hop_data.labels_col));

regex_string = '%*s scores for Exp. %d_%f_HOP_%s %s';
for i = 1 : length(hop_data.labels_col)
    C = textscan(hop_data.labels_col{i}, regex_string);
    hop_data.cmd(i) = C{1};
    hop_data.conc(i) = C{2};
    hop_data.group(i) = C{3};
    if ~isempty(C{4})
        hop_data.zscore(i) = C{4};
    end
end

% Only retain the z-scores
inds = find(~strcmp('z-score', hop_data.zscore));
hop_data.labels_col(inds) = [];
hop_data.cmd(inds) = [];
hop_data.conc(inds) = [];
hop_data.group(inds) = [];
hop_data.zscore(inds) = [];
hop_data.data(:,inds) = [];


% Map compound IDs to common names (in the few cases where that's possible)
[data.txt, data.num, data.raw] = xlsread('/Users/Anastasia/Laboratory/Datasets/Phenotypes/2014_Hoepfner~Movva/Table_S1.xlsx','Reference Substances known MoA');
cmp.id = cell2mat(data.raw(2:end,1));
cmp.name = data.raw(2:end,2);
inds = find(isnan(cmp.id));
cmp.id(inds) = [];
cmp.name(inds) = [];

[Lia, Locb] = ismember(hop_data.cmd, cmp.id);
hop_data.cmd_name = cell(size(hop_data.cmd));
hop_data.cmd_name(Lia) = cmp.name(Locb(Lia));

% Check the ORFs
inds = cellfun(@isnumeric, hop_data.labels_row);
hop_data.labels_row(inds) = [];

hop_data.labels_row = cellfun(@(x) regexprep(x,'"',''), hop_data.labels_row, 'UniformOutput',0);    % remove the unnecessary quotes from strings
hop_data.labels_row = cellfun(@strtrim, hop_data.labels_row,'UniformOutput',0);

hop_data.labels_row = upper(hop_data.labels_row);

regexp_orf = 'Y[A-P][RL][0-9]{3}[CW](-[ABC])*';
inds = find(cellfun(@isempty, regexp(hop_data.labels_row, regexp_orf)));
hop_data.labels_row(inds) = [];
hop_data.data(inds,:) = [];

% Print sample for SAFE analysis
inds = find(~cellfun(@isempty, hop_data.cmd_name));
ucmd = unique(hop_data.cmd_name(inds));
[~,ind1,ind2] = intersect(hop_data.cmd_name(inds), ucmd);

fid = fopen('~/Laboratory/Datasets/Phenotypes/2014_Hoepfner~Movva/hoepfner_movva_2014_hop_known.txt','w');
print_matrix_file(fid, hop_data.labels_row, hop_data.cmd_name(inds(ind1)), hop_data.data(:,inds(ind1)));
fclose(fid);




