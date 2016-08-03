classdef sParameterPlot
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fh
        plotType
    end
    
    methods
        % constructor
        function self = sParameterPlot(fh, plotType)
            if nargin < 1 || fh < 1
                self.fh = figure();
            else
                self.fh = figure(fh);
            end
            if nargin < 2
                self.plotType = 'dB';
            elseif strcmpi(plotType,'db')
                self.plotType = 'dB';
            elseif ~isempty(regexpi(plotType,'^mag')) || strcmpi(plotType,'abs')
                self.plotType = 'abs';
            elseif ~strcmpi(plotType,'phase') || ~strcmpi(plotType,'angle')
                self.plotType = 'angle';
            end
        end
    end
    
end

