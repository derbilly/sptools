function [NetData,Frequency] = get_citi(spfile)
% GET_CITI imports a CITI file into
% a frequency vector and a complex three-dimensional scattering
% parameter matrix.
%
% William R. Peters 12/9/2004
%
% USAGE: [NetData,Frequency] = get_citi(spfile)
%  where
%   Frequency(k)                kth frequency
%   NetData(n,m,k)              Snm at kth frequency (complex)
%   spfile                      CITI file

% Open file for read ====================
FID = fopen(spfile);
SourceFile=spfile;

% Get CITIFILE version info ==================
nextline = GetNextLine(FID);
[flag,version] = strtok(nextline);
if ~strcmp(upper(flag),'CITIFILE') % check for correct syntax
    error('unexpected field, CITIFILE expected');
end

% Get name ==================
nextline = GetNextLine(FID);
[flag,name] = strtok(nextline);
if ~strcmp(upper(flag),'NAME') % check for correct syntax
    error('unexpected field, NAME expected');
end

% Get independent variable info ==================
nextline = GetNextLine(FID);
[flag,remainder] = strtok(upper(nextline));
if ~strcmp(flag,'VAR') % check for correct syntax
    error('unexpected field, VAR expected');
end
[varname,remainder] = strtok(remainder);
[varform,remainder] = strtok(remainder);
varsteps = str2num(strtok(remainder));
if ~strcmp(varname,'FREQ') | ~strcmp(varform,'MAG') | ~isnumeric(varsteps)   
    error('VAR format error')
end

% Get data field info
nd=0; md=0; datakey=[]; dataform=[];
nextline = GetNextLine(FID);
[flag,remainder] = strtok(upper(nextline));
while ~strcmp(flag,'SEG_LIST_BEGIN')
    if strcmp(flag,'DATA') & strcmp(remainder(1:3),' S[')
        nd=nd+1;
        [dkey,dform] = strtok(remainder);
        datakey{nd}=dkey;
        dataform{nd}=dform;
    elseif strcmp(flag,'DATA')
        md=md+1;
        % add capability here
    else
        error('unexpected field, expected DATA');
    end
    nextline = GetNextLine(FID);
    [flag,remainder] = strtok(upper(nextline));
end
if rem(sqrt(nd),1)
    error('error in reading data format')
end
N=sqrt(nd);
% check data format lines
for p=1:N
    for q=1:N
        if ~strcmp(datakey{(p-1)*N+q},sprintf('S[%d,%d]',p,q))
            error(sprintf('Data key error %d,%d',p,q))
        end
    end
end
% check data format consistency
if ~strcmp(dataform,dataform{1})
    error('data format inconsistency')
end

% get frequency data
nextline = GetNextLine(FID);
[flag,remainder] = strtok(upper(nextline));
if ~strcmp(flag,'SEG') % check for correct syntax
    error('unexpected field, expected SEG');
end
freqstep = sscanf(remainder,'%f %f %f');
if freqstep(3)~=varsteps
    error('VAR SEG inconsistency error')
end

Frequency = linspace(freqstep(1),freqstep(2),freqstep(3));
nextline = GetNextLine(FID);
[flag,remainder] = strtok(upper(nextline));
if ~strcmp(flag,'SEG_LIST_END') % check for correct syntax
    error('unexpected field, expected SEG_LIST_END');
end

NetData=zeros(N,N,freqstep(3));
for p=1:N
    for q=1:N
        nextline = GetNextLine(FID);
        [flag,remainder] = strtok(upper(nextline));
        if ~strcmp(flag,'BEGIN') % check for correct syntax
            error('unexpected field, expected BEGIN');
        end
        for r=1:freqstep(3)
           nextline = GetNextLine(FID);
           datc = sscanf(nextline,'%f,%f');
           NetData(p,q,r) = datc(1) + i*datc(2);
        end
        nextline = GetNextLine(FID);
        [flag,remainder] = strtok(upper(nextline));
        if ~strcmp(flag,'END') % check for correct syntax
            error('unexpected field, expected END');
        end
    end
end    
fclose(FID);
    
function out = GetNextLine(FID)
% get the next significant line
goodline = 0;
while ~goodline
    out = fgetl(FID);
    if ~ischar(out)
        if out==-1
            return
        end
    end
    if IsComment(out)
        goodline=0;
    else
        goodline=1;
    end
end


function out = IsComment(line)
% determine if line is a comment or blank
testline = strrep(line,' ','');
if isempty(testline)
    out=1;
elseif testline(1)=='#'
    out=1;
elseif strcmp(upper(strtok(line)),'COMMENT')
    out=1;
else
    out=0;
    line
end