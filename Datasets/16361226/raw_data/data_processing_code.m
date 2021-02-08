%% Cai~Becker, 2006
function FILENAMES = code()

addpath(genpath('../../Yeast-Matlab-Utils/'));

FILENAMES = {};
cai_becker_2006.pmid = 16361226;

% MANUAL. Download the list of dataset ids and standard names from
% the paper's page on www.yeastphenome.org & save the file to ./extras

% Load the list
[FILENAMES{end+1}, d] = read_data('textread', ['./extras/YeastPhenome_' num2str(cai_becker_2006.pmid) '_datasets_list.txt'],'%d %s','delimiter','\t');
datasets.id = d{1};
datasets.standard_name = d{2};

%% Load the data

flds = dir('./raw_data/data_processing/Sceeningdeletionmutant_rawdata/');
flds(1:3) = [];

% PLATE MAPS: Get the data from the Excel files
map = [];
for i = 1 : length(flds)
    tmp = regexp(flds(i).name, '_', 'split');
    plate = str2num(tmp{1});
    
    files = dir(['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
    f = find(~cellfun(@isempty, strfind({files.name},'all.xls')));
    filename = ['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f).name];
    [~,sheets] = xlsfinfo(filename);
    
    for sht = 1 : 3
        
        [FILENAMES{end+1}, mapdata.raw] = read_data('xlsread',filename, sheets{sht});
        
        if length(find(cellfun(@isnumeric,mapdata.raw(:,4)))) < 10
            map = [map; mapdata.raw(:,1:5)];
            break;
        end
    end
end

inds = find(~cellfun(@ischar, map(:,1)));
map(inds,1) = {'EMPTY'};

map(:,1) = upper(map(:,1));

% PLATE MAPS: A few rows to eliminate
original = {'STRAIN','FILTER'};
for i = 1 : length(original)
    inds = find(strncmp(original{i}, map(:,1), length(original{i})));
    map(inds,:) = [];
end
inds = find(cellfun(@isnumeric, map(:,4))); % Rows with "NaN" coordinates
map(inds,:) = [];
inds = find(strcmp(' ', map(:,4))); % Rows with " " (empty) coordinates
map(inds,:) = [];


% PLATE MAPS: A few rows to rename
original = {'NONE','PTR1','PTR2','CUP9'};
new = {'EMPTY','YGR184C','YKR093W','YPL177C'};
for i = 1 : length(original)
    inds = find(strncmp(original{i}, map(:,1), length(original{i})));
    map(inds,1) = new(i);
end

% PLATE MAPS: Change row letters into numbers
letters = 'ABCDEFGH';
for i = 1 : size(map, 1)
    map{i,4} = strfind(letters, map{i,4});
end

% PLATE MAPS: Keep ORFs and background as cell array, while plate, row &
% columns as numbers
map_txt = map(:,1:2);
map_num = cell2mat(map(:,3:5));


% % TIMEPOINTS (1): TO DO ONLY ONCE. Since this information is not standardized, it is necessary to generate a SUMMARY FILE & process it manually
% fid = fopen('./raw_data/data_processing/summary.txt','w');
% for i = 1 : length(flds)
%     
%     tmp = regexp(flds(i).name, '_', 'split');
%     plate = str2num(tmp{1});
%     if length(replicates) < plate
%         replicates(plate) = 1;
%     else
%         replicates(plate) = replicates(plates) + 1;
%     end
%     
%     files = dir(['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
%     
%     f = find(cellfun(@isempty, strfind({files.name},'.dat'))==0);
%     
%     for j = 1 : length(f)
%         
%         filename = ['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f(j)).name];
%         
%         [FILENAMES{end+1}, C] = dataread('textread',filename, '%s', 'delimiter', '\n');
%         
%         if length(C) >= 9
%             header = C(1);
%             fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
%         end
%         
%         if length(C) >= 19
%             header = C(11);
%             fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
%         end
%         
%         if length(C) >= 29
%             header = C(21);
%             fprintf(fid,'%s\t%s\t%d\t%d\t%s\n', flds(i).name, files(f(j)).name, plate, replicates(plate_red), header{h});
%         end
%         
%         clear C;
%     end
%     
% end
% fclose(fid);

% TIMEPOINTS (2): The summary.txt file was processed manually to extract information about 1) Medium, 2) Date, 3) Timepoint.

% TIMEPOINTS (3): Load in the new summary2.txt file
[FILENAMES{end+1}, C] = read_data('importdata','./raw_data/data_processing/summary2.txt', '\t');
fields = regexp(C{1}, '\t', 'split');
fields(end) = {'TMP'};

for i = 2 : length(C)
    tmp = regexp(C{i}, '\t', 'split');
    for j = 1 : length(fields)
        summary.(fields{j})(i,1) = tmp(j);
    end
end

% DATA (1): Go over all the data and save it properly

flds = dir('./raw_data/data_processing/Sceeningdeletionmutant_rawdata/');
flds(1:3) = [];

data = struct('plateid',[0], 'medium',{'x'},'plate',[0],'date',{'x'},'hour',{'x'},'timepoint',[0],'replicates',[0],'orfs',{'x'},'rows',[0],'cols',[0],'data',[0],'folder',{'x'}, 'file',{'x'},'header',{'x'});

plateid = 0;    %unique plate identifier

for i = 1 : length(flds)
    
    % To identify a set of data you need folder, file and header.
    
    inds1 = find(strcmp(flds(i).name, summary.Folder)); % FOLDER
    
    tmp = regexp(flds(i).name, '_', 'split');
    plate = str2num(tmp{1});
    
    % Find the ORFs
    inds = find(map_num(:,1) == plate);
    orfs = map_txt(inds,1);
    rows = map_num(inds,2);
    cols = map_num(inds,3);
    
    inds_data = sub2ind([8 12], rows, cols);
    
    files = dir(['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/']);
    
    f = find(cellfun(@isempty, strfind({files.name},'.dat'))==0);
    
    for j = 1 : length(f)
        
        inds2 = find(strcmp(files(f(j)).name, summary.File)); % FILE
        
        filename = ['./raw_data/data_processing/Sceeningdeletionmutant_rawdata/' flds(i).name '/' files(f(j)).name];
        
        [FILENAMES{end+1}, C] = read_data('textread',filename, '%s', 'delimiter', '\n');
        
        data_row_num = [1 11 21];
        
        for d = 1 : 3
            
            if length(C) > data_row_num(d)
                
                header = C(data_row_num(d));
                inds3 = find(strcmp(header{1}, summary.Header)); % HEADER
                inds = intersect(intersect(inds1,inds2),inds3); % IDENTIFY DATA
                
                [FILENAMES{end+1}, A] = read_data('dlmread',filename, '\t', [data_row_num(d) 0 data_row_num(d)+7 11]);
                
                if ~isempty(summary.Date{inds}) & ~strcmp('#VALUE!', summary.Hour{inds})
                    
                    plateid = plateid+1;
                    
                    data.plateid = [data.plateid; zeros(length(orfs),1)+plateid];
                    data.medium = [data.medium; repmat(summary.Medium(inds),length(orfs),1)];
                    data.plate = [data.plate; repmat(plate, length(orfs),1)];
                    data.date = [data.date; repmat(summary.Date(inds), length(orfs),1)];
                    data.hour = [data.hour; repmat(summary.Hour(inds), length(orfs),1)];
                    data.timepoint = [data.timepoint; repmat(summary.Timepoint(inds), length(orfs),1)];
                    data.replicates = [data.replicates; repmat(summary.Replicates(inds), length(orfs),1)];
                    data.orfs = [data.orfs; orfs];
                    data.rows = [data.rows; rows];
                    data.cols = [data.cols; cols];
                    data.data = [data.data; A(inds_data)];
                    data.folder = [data.folder; repmat(summary.Folder(inds), length(orfs),1)];
                    data.file = [data.file; repmat(summary.File(inds), length(orfs),1)];
                    data.header = [data.header; repmat(header, length(orfs),1)];
                    
                end
            end
        end
    end
end


fields = fieldnames(data);
for i = 1 : length(fields)
    data.(fields{i})(1) = [];   % delete the first initiating value.
end

% DATA (2): Transform a few values into numbers
data.timepoint = cellfun(@str2num, data.timepoint); % timepoints
data.replicates = cellfun(@str2num, data.replicates); % replicates
[tmp,ia,data.medium_num] = unique(data.medium); % medium

[tmp_u,ia,ib] = unique(strcat(data.date, {' '}, data.hour)); % date & time
tmp_u_vc = datevec(tmp_u);
tmp_u_vc(tmp_u_vc(:,2)>1,1) = 2003; tmp_u_vc(tmp_u_vc(:,2)==1,1) == 2004;
tmp_u_nm = datenum(tmp_u_vc);
data.datenum = tmp_u_nm(ib);

% DATA (3): Calculate proper time intervals between timepoints
plate_rep_medium = unique([data.plate data.replicates data.medium_num],'rows');
for i = 1 : size(plate_rep_medium,1)
    
    inds = find(data.plate == plate_rep_medium(i,1) & data.replicates == plate_rep_medium(i,2) & data.medium_num == plate_rep_medium(i,3));
    [timepoints,ia,ib] = unique(data.timepoint(inds));
    datetime = data.datenum(inds(ia));
    timeintervals = [0; diff(datetime)] *24;    % number of hours between consecutive timepoints
    
    data.timeintervals(inds) = timeintervals(ib);
    
end


% % QC: Visualize the data to see if there're any biases.
% inds = find(data.plateid == randsample(unique(data.plateid),1)); % grab a random plate
% otherplates = find(data.plate == data.plate(inds(1)) & data.replicates == data.replicates(inds(1)) & strcmp(data.medium{inds(1)}, data.medium)); % find the other timepoints for this plate, this medium and this replicate
%
% plateids = unique(data.plateid(otherplates));
% figure();
% for i = 1 : length(plateids)
%     subplot(3,3,i);
%     plate_map = zeros(8,12);
%     inds = find(data.plateid == plateids(i));
%     data_inds = sub2ind([8 12], data.rows(inds), data.cols(inds));
%     plate_map(data_inds) = data.data(inds);
%
%     imagesc(plate_map,[0 0.3]);
%     title([data.medium{inds(1)} ', Plate ' num2str(data.plate(inds(1))) ', Replicate ' num2str(data.replicates(inds(1))) ', Date/time ' data.date{inds(1)} '/' data.hour{inds(1)}]);
%     colorbar;
%     axis image;
% end

% DATA (4): Combine all timepoints into a single measure. Arbitrary = Area
% under the growth curve

% Prepare summary matrix ORFS x Medium x Replicates
all_orfs = unique(data.orfs);
all_media = unique(data.medium);

plate_rep_medium = unique([data.plate data.replicates data.medium_num],'rows');
data_matrix = zeros(length(all_orfs), length(all_media),size(plate_rep_medium,1))+NaN;
curr_rep = 1;
for i = 1 : size(plate_rep_medium,1)
    inds = find(data.plate == plate_rep_medium(i,1) & data.replicates == plate_rep_medium(i,2) & data.medium_num == plate_rep_medium(i,3));
    [timeintervals,ia,ib] = unique(data.timeintervals(inds));
    [orfs,i1,i2] = unique(data.orfs(inds));
    auc = zeros(length(orfs),1);
    for o = 1 : length(orfs)
        inds2 = find(strcmp(orfs{o}, data.orfs(inds)));
        [t,ix] = sort(data.timepoint(inds(inds2)));
        dt = data.data(inds(inds2(ix)))-data.data(inds(inds2(ix(1))));  % difference in growth between now and the first timepoint
        ds = data.timeintervals(inds(inds2(ix)))';
        auc(o) = nansum(ds.*dt);
        if isnan(auc(o))
            break;
        end
    end
    % Normalize by the average of the entire plate
    auc_norm = auc ./ nanmean(auc);
    
    [t,ind1,ind2] = intersect(all_orfs, orfs);
    data_matrix(ind1,plate_rep_medium(i,3),curr_rep) = auc_norm(ind2);
    curr_rep = curr_rep + 1;
end

% Check how many replicates ORFs have in general
n = sum(~isnan(data_matrix),3);

% Average data for ORFs across all replicates
data_matrix2 = nanmean(data_matrix,3);

%%%% -> BEGIN: printed to TXT and continued in python from this point on

fid = fopen('processed_data.txt','w');
col_headers = {'108','245','246'}';
write_matrix_file(fid, all_orfs, col_headers, data_matrix2);
fclose(fid);

%%%% -> END

% % Checks on ORFs
% all_orfs = clean_orf(all_orfs);
% 
% inds = find(~is_orf(all_orfs));
% disp(all_orfs(inds));
% 
% all_orfs(inds) = [];
% data_matrix2(inds,:) = [];
% 
% % MANUAL. Get the dataset ids corresponding to each dataset (in order)
% % Multiple datasets (e.g., replicates) may get the same id, which can then
% % be used to average them out
% hit_data_ids = [108 245 246]';
% 
% %% Prepare final dataset
% 
% % Match the dataset ids with the dataset standard names
% [~,ind1,ind2] = intersect(datasets.id, hit_data_ids);
% hit_data_names = cell(size(hit_data_ids));
% hit_data_names(ind2) = datasets.standard_name(ind1);
% 
% cai_becker_2006.orfs = all_orfs;
% cai_becker_2006.ph = hit_data_names;
% cai_becker_2006.data = data_matrix2;
% cai_becker_2006.dataset_ids = hit_data_ids;
% 
% %% Save
% 
% save('./cai_becker_2006.mat','cai_becker_2006');
% 
% %% Print out
% 
% fid = fopen('./cai_becker_2006.txt','w');
% write_matrix_file(fid, cai_becker_2006.orfs, cai_becker_2006.ph, cai_becker_2006.data);
% fclose(fid);
% 
% %% Save to DB (admin)
% 
% addpath(genpath('../../Private-Utils/'));
% if exist('save_data_to_db.m')
%     res = save_data_to_db(cai_becker_2006)
% end

end