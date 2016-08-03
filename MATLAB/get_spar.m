function [NetData,Frequency] = get_spar(spfile)
% GET_SPAR imports a CITI or Touchstone file into
% a frequency vector and a complex three-dimensional scattering
% parameter matrix.
%
% William R. Peters 12/13/2004
%
% USAGE: [NetData,Frequency] = get_spar(spfile)
%  where
%   Frequency(k)           kth frequency
%   NetData(n,m,k)               Snm at kth frequency (complex)
%   spfile                      CITI or Touchstone file

% Open file for read ====================
FID = fopen(spfile);
% Determine S parameter format ==================
nextline = GetNextLine(FID);
[flag,version] = strtok(nextline);
if strcmp(upper(flag),'CITIFILE')
    fclose(FID);
    [NetData,Frequency] = get_citi(spfile);
elseif strcmp(upper(flag),'#')
    fclose(FID);
    [NetData,Frequency] = get_snp(spfile);
else
    fclose(FID);
    error('unrecognized format');
end    
    
    
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
line = strrep(line,' ','');
if isempty(line)
    out=1;
elseif line(1)=='!'
    out=1;
else
    out=0;
end