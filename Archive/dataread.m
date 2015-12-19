function [filenames, DATA] = dataread(varargin) %type, filename, ...
%this function replaces all the 
filenames = varargin{2};
type = varargin{1};

switch type
    case 'textread' %be able to grab any number of outputs
        fid = fopen(varargin{2},'r');
        DATA = textscan(fid, varargin{3:end});
        if numel(DATA) == 1
            DATA = DATA{1};
        end
        fclose(fid);
    case 'textscan'
        fid = fopen(varargin{2},'r');
        DATA = textscan(fid, varargin{3:end});
        if numel(DATA) == 1
            DATA = DATA{1};
        end
        fclose(fid);
    case 'xlsread' 
        [~, ~, DATA] = xlsread(varargin{2:end}); %Do I only accept the final output?
    case 'dlmread'
        DATA = dlmread(varargin{2:end}); %This should work fine.
    case 'read_matrix_file'
        DATA = read_matrix_file(varargin{2:end}); %This works perfectly fine
    case 'importdata'
        DATA = importdata(varargin{2:end}); %This should work fine.
    otherwise
        error('Not Valid Function');
end



end