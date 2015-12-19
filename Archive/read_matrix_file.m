function D = read_matrix_file(input_file, numrowheaders, numcolheaders)

%
% [labels_row, labels_col, data] = read_matrix_file(input_file, num_row_headers, num_col_headers)
% TREEVIEW FILE: num_row_headers = 2, num_col_headers = 2
% 

if nargin < 3 || isempty(numrowheaders)
    numrowheaders = 1;
end

if nargin < 3 || isempty(numcolheaders)
    numcolheaders = 1;
end

D.numrowheaders = numrowheaders;
D.numcolheaders = numcolheaders;

% D.data = dlmread(input_file,'\t', numrowheaders, numcolheaders);

fid = fopen(input_file,'r');
if fid == -1
    fprintf('Unable to open file %s\n', input_file);
    return;
end


% Get column headers
D.labels_col = {};
for i = 1 : numcolheaders
    line = fgetl(fid);
    tmp = textscan(line,'%s','delimiter','\t');
    colheaders = tmp{1}(numrowheaders+1:end);
    if length(colheaders) > size(D.labels_col,1)
        tmp = D.labels_col;
        D.labels_col = cell(length(colheaders),size(D.labels_col,2));
        D.labels_col(1:size(tmp,1),1:size(tmp,2)) = tmp;
        D.labels_col(:,i) = colheaders;
    else
        D.labels_col(1:length(colheaders),i) = colheaders;
    end
end

% Get row headers and data
fmt = [repmat('%s ', 1, numrowheaders) repmat('%f ', 1, length(D.labels_col))];

D.labels_row = cell(1,numrowheaders);
D.data = [];

r = 1;
while 1
    line = fgetl(fid);
    if line == -1 
        break; 
    end
    
    if mod(r,100) == 0
        fprintf('Read %d lines\n', r);
    end
       
    t = textscan(line,fmt,'delimiter','\t');
    t(find(cellfun(@isempty,t))) = {NaN};
    
    D.labels_row(r,:) = cat(2,t{1:numrowheaders});
    D.data(r,:) = cell2mat(t(numrowheaders+1:end));
    
    r = r + 1;
    
end