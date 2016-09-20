classdef SParameters < FrequencyDomainData
    %SParameters - Base class for S-parameters
    %   Detailed explanation goes here
    
    properties
        label
    end
    
    properties (SetAccess = protected)
        numPorts
        frequency
        portZ
        S
        dataFile
    end
    
    
    methods
        % constructor
        function self = SParameters(touchstoneFile)
            self@FrequencyDomainData();
            if ~exist('touchstoneFile')
                [filename,path] = uigetfile('*.s*p');
                touchstoneFile = [path,filename];
            end
            [self.S,self.frequency,self.numPorts,self.portZ] = get_snp(touchstoneFile);
            self.dataFile = touchstoneFile;
            self.label = regexprep(self.dataFile,'\.s[0-9]*p$','');
            self.label = regexprep(self.label,'.*\\','');
            self.label = regexprep(self.label,'_',' ');
        end
        
        % identify
        function identify(self)
            disp(sprintf('%i port S-Parameter object defined from %g Hz to %g Hz.',self.numPorts,self.frequency(1),self.frequency(end)));
        end
        
       
        % more methods
        function self = reorderPorts(self,portOrder)
            % check inputs!!!!
            self.S = self.S(portOrder,portOrder,:);
        end
        
        function self = resampleFrequency(self,frequencyNew)
            % only interpolation so far
            Snew = zeros([self.numPorts,self.numPorts,length(frequencyNew)]);
            for n = 1:self.numPorts
                for m = 1:self.numPorts
                    magnitude0 = squeeze(abs(self.S(n,m,:)));
                    magnitude = interp1(self.frequency,magnitude0,frequencyNew,'pchip');
                    phase0 = squeeze(angle(self.S(n,m,:)));
                    phase = interp1(self.frequency,unwrap(phase0),frequencyNew,'pchip');
                    Snew(n,m,:) = magnitude.*exp(1i.*phase);
                end
            end
            self.frequency = frequencyNew;
            self.S = Snew;
        end
        
        function plotAll(self,matrix)
            figure()
            if ~exist('matrix')
                matrix = 'S';
            end
            N = size(self.(matrix),1);
            for n = 1:N
                for m = 1:N
                    subplot(N,N,(n-1)*N+m)
                    plot(self.frequency.*1e-9,dB(squeeze(self.(matrix)(n,m,:))))
                    xlabel('GHz')
                    ylabel('dB')
                    title(sprintf('%s(%d,%d)',matrix,n,m))
                    grid on
                end
            end
        end
        
        % touchstone exporter
        
        % getS elaborate this
        function out = getS(self,n,m)
            out  = squeeze(self.S(n,m,:));
        end
        
    end
    
end

%% Local functions
%% Touchstone Importer
function [NetData,frequency,numPorts,portZ] = get_snp(spfile)
% GETSNP imports a touchstone file into
% a frequency vector and a complex three-dimensional scattering
% parameter matrix.
%
% William R. Peters 7/26/04, 2/18/05, 6/2/16
%
% USAGE: [NetData,frequency,numPorts,portZ] = get_snp(spfile)
%  where
%   frequency(k)    kth frequency
%   NetData(n,m,k)  Snm at kth frequency (complex)
%   numPorts        number of ports
%   portZ           port impedance
%   spfile          touchstone file
%

% Open file for read ====================
FID = fopen(spfile);
SourceFile=spfile;
% Get format line info ==================
nextline = GetNextLine(FID);
[formsign,remainder] = strtok(upper(nextline));
if formsign~='#' % check for correct syntax
    error('error in sp file format');
end
[FreqUnit,remainder] = strtok(remainder); % get frequency units
[NetForm,remainder] = strtok(remainder); % check that s params are specified
if NetForm~='S'
    error('need s parameter format');
end
[ComplexForm,remainder] = strtok(remainder); % get complex form (MA,RI,DB)
[R,remainder] = strtok(remainder); % check for correct syntax
if R~=R
    error('format line error (R not found)');
end
[portZ,remainder] = strtok(remainder); % check that port impedance is 50 ohms
portZ = str2num(portZ);
%if str2num(portZ) ~= 50
%    error('port impedance must be 50 ohms')
%elseif length(remainder) ~= 0
%    error('must have real port impedance')
%end
% get data
n=0;
nextline = GetNextLine(FID);
while nextline ~= -1
    n=n+1;
    CellData{n} = nextline;
    nextline = GetNextLine(FID);
end
% find number of lines per frequency point
Nlines=1;
if rem(length( sscanf(CellData{Nlines},'%f') ) , 2) == 0
    error('even number of entries in 1st line')
end
while rem(length( sscanf(CellData{Nlines+1},'%f') ) , 2) == 0
    Nlines = Nlines + 1;
end
if rem( length(CellData), Nlines ) ~= 0
    error('incorrect number of lines in file')
end
% create data matrix with one line per frequency point
for n = 1:length(CellData)/Nlines
    DataMatrix0 = sscanf(CellData{(n-1)*Nlines+1},'%f');
    for m = 1:Nlines-1
        DataMatrix0 = [DataMatrix0 ; sscanf(CellData{(n-1)*Nlines+1+m},'%f')];
    end
    DataMatrix(n,:) = DataMatrix0';
end
% convert data matrix into complex 3-D S matrix and frequency
% (from old function)
%
% frequency multiplier
if strcmp(FreqUnit,'HZ')
    FreqMult=1;
elseif strcmp(FreqUnit,'KHZ')
    FreqMult=1e3;
elseif strcmp(FreqUnit,'MHZ')
    FreqMult=1e6;
elseif strcmp(FreqUnit,'GHZ')
    FreqMult=1e9;
else
    error('frequency unit problem')
end
% create frequency vector
frequency = DataMatrix(:,1)*FreqMult;
% create s-parameter matrix
Nf = length(frequency);
numPorts = sqrt((size(DataMatrix,2) - 1)/2);
NetData = zeros(numPorts,numPorts,Nf);
nn = 0;
if ComplexForm=='RI' % real imaginary
    for rr = 1:numPorts
        for cc = 1:numPorts
            nn = nn + 1;
            NetData(rr,cc,:) = DataMatrix(:,nn*2) + i*DataMatrix(:,nn*2+1);
        end
    end
elseif ComplexForm=='MA' % magnitude angle
    for rr = 1:numPorts
        for cc = 1:numPorts
            nn = nn+1;
            NetData(rr,cc,:) = DataMatrix(:,nn*2).*exp(i*(2*pi/360)*DataMatrix(:,nn*2+1));
        end
    end
elseif ComplexForm=='DB' % dB angle
    for rr = 1:numPorts
        for cc = 1:numPorts
            nn = nn+1;
            NetData(rr,cc,:) = 10.^(DataMatrix(:,nn*2)/20).*exp(i*(2*pi/360)*DataMatrix(:,nn*2+1));
        end
    end
end
fclose(FID);
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
end


