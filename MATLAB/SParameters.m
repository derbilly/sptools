classdef SParameters
    %SPARAM Summary of this class goes here
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
            [self.S,self.frequency,self.numPorts,self.portZ] = get_snp(touchstoneFile);
            self.dataFile = touchstoneFile;
            self.label = regexprep(regexprep(self.dataFile,'\.s[0-9]*p$',''),'_',' ');
        end
        
        % identify
        function identify(self)
            disp(sprintf('%i port S-Parameter object defined from %g Hz to %g Hz.',self.numPorts,self.frequency(1),self.frequency(end)));
        end
        
        %         % setters
        %         %function self = set.S(self,value)
        %         %    self.S = value;
        %         %end
        
        % getters
        %         function out = get.frequency(self)
        %             out = self.frequency;
        %         end
        %
        %         function out = get.S(self)
        %             out = self.S;
        %         end
        %
        %         function out = get.portZ(self)
        %             out = self.portZ;
        %         end
        %
        %         function out = get.label(self)
        %             out = self.label;
        %         end
        %
        %         function out = get.numPorts(self)
        %             out = self.numPorts;
        %         end
        
        % more methods
        function powersum = checkPassivity(self)
            % checks for passivity
            powersum = zeros(length(self.frequency),self.numPorts);
            for porti = 1:self.numPorts
                for portj = 1:self.numPorts
                    powersum(:,portj) = powersum(:,portj) + squeeze(self.S(porti,portj,:).*conj(self.S(porti,portj,:)));
                end
            end
        end
        
        
        function addToPlot(self,varargin)%figure_num,item,index,plotSpec
            %defaults
            if length(varargin) < 4
                plotSpec = '-';
            else
                plotSpec = varargin{4};
            end
            if length(varargin) < 3
                index = [1,1];
            else
                index = varargin{3};
            end
            if length(varargin) < 2
                item = 'S';
            else
                item = varargin{2};
            end
            if length(varargin) < 1 || varargin{1} == 0
                0;
            else
                figure(varargin{1});
            end
            [~,~,~,OUTM] = legend;
            hold on;
            plot(self.frequency*1e-9,squeeze(dB(self.(item)(index(1),index(2),:))),plotSpec);
            hold off;
            OUTM{end+1} = [ self.label ' ' item '_{' num2str(index(1)) num2str(index(2)) '}' ];
            legend(OUTM,'Location','Best');
        end
        
        
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
                    unwrappedPhase0 = unwrapPhase(phase0);
                    phase = interp1(self.frequency,unwrappedPhase0,frequencyNew,'pchip');
                    Snew(n,m,:) = magnitude.*exp(1i.*phase);
                end
            end
            self.frequency = frequencyNew;
            self.S = Snew;
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

function angleOut = unwrapPhase(angleIn)

angleOut = zeros(size(angleIn));
angleOut(1)=angleIn(1);
m=0;
for n = 2:length(angleIn)
    if angleIn(n)-angleIn(n-1) > pi
        m=m+1;
    elseif angleIn(n)-angleIn(n-1) < -pi
        m=m-1;
    end
    angleOut(n)=angleIn(n)-m*2*pi;
end
end
