function legendString = plotSparam( plotSpec, plotType, titleString )
%PLOTSPARAM Summary of this function goes here
%   Detailed explanation goes here


if nargin==1
    plotType = 'dB';
elseif strcmpi(plotType,'db')
    plotType = 'dB';
elseif ~isempty(regexpi(plotType,'^mag')) || strcmpi(plotType,'abs')
    plotType = 'abs';
elseif ~strcmpi(plotType,'phase') || ~strcmpi(plotType,'angle')
    plotType = 'angle';
end

plotString = {};
legendString = {};
plotItem = 0;
for n = 1:size(plotSpec,1)
    dataName = plotSpec{n,2};
    indices = plotSpec{n,3};
    for m=1:size(indices,1)
        plotItem = plotItem+1;
        plotString{plotItem} = [ 'plotSpec{' num2str(n) ',1}.frequency*1e-9,squeeze(' plotType '(plotSpec{' num2str(n) ',1}.' dataName '(' num2str(indices(m,1)) ',' num2str(indices(m,2)) ',:)))' ];
        legendString{plotItem} = [ plotSpec{n,1}.label ' ' dataName  num2str(indices(m,1)) num2str(indices(m,2))  ];
    end
end

% eval(['plot(' strjoin(plotString,',') ')' ]);
% rewrite for MATLAB 2012 or earlier without strjoin
for n = 1:length(plotString)-1
    plotString{n} = [ plotString{n},',' ];
end
plotString  = cell2mat(plotString);
%plotString = regexprep(plotString,',$',''); %remove last comma
eval(['plot(' plotString ')' ]);

if nargin >= 3
    title(titleString);
end
xlabel('GHz')
if strcmp(plotType,'dB')
    ylabel('dB')
elseif strcmp(plotType,'abs')
    ylabel('magnitude')
elseif strcmp(plotType,'angle')
    ylabel('degrees')
end

legend(legendString)
grid on
end

function out = dB(x)
% Computes decibels without signal processing toolbox
% = 20*log10(abs(x));
out = 20*log10(abs(x));
end

