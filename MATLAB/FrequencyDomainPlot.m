classdef FrequencyDomainPlot < DataPlot
    properties
        dataFormat
        title
        dataList
        frequencyList
        labelList
        plotSpecList
    end
    
    methods
        function self = FrequencyDomainPlot(dataFormat,title)
            if exist('dataFormat','var')
                self.dataFormat = dataFormat;
            else
                self.dataFormat = 'dB';
            end
            if exist('title','var')
                self.title = title;
            else
                self.title = '';
            end
            self.dataList = {};
            self.frequencyList = {};
            self.labelList = {};
            self.plotSpecList = {};
        end
        
        function self = addItem(self,data,arg2,arg3,lineStyle)
            if ~exist('arg2','var')
                arg2 = 'SDD';
            end
            if ~exist('arg3','var')
                arg3 = [1,1];
            end
            if isa(data,'FrequencyDomainData')
                if isa(data,'SParameters')
                    self.dataList{end+1} = squeeze(data.(arg2)(arg3(1),arg3(2),:));
                    self.frequencyList{end+1} = data.frequency;
                    self.labelList{end+1} = sprintf('%s %s(%d,%d)',data.label,arg2,arg3(1),arg3(2)); % fix
                    if ~exist('lineStyle','var')
                        lineStyle = '';
                    end
                    self.plotSpecList{end+1} = lineStyle;
                elseif isa(data,'SpecLine')
                    self.dataList{end+1} = data.specLine;
                    self.frequencyList{end+1} = data.frequency;
                    if ~exist('lineStyle','var')
                        lineStyle = '--';
                    end
                    self.plotSpecList{end+1} = lineStyle;
                    self.labelList{end+1} = data.specItem; % fix
                end
            else
                self.dataList{end+1} = data;
                self.frequencyList{end+1} = arg2;
                self.labelList{end+1} = label;
                if ~exist('lineStyle','var')
                    lineStyle = '';
                end
                self.plotSpecList{end+1} = lineStyle;
            end
        end
        
        function generatePlot(self)
            figure();
            hold on
            % add limits
            for n = 1:length(self.dataList)
                if strcmp(self.dataFormat,'dB')
                    plot(self.frequencyList{n}/1e9,dB(self.dataList{n}),self.plotSpecList{n});
                elseif strcmp(self.dataFormat,'radians')
                    plot(self.frequencyList{n}/1e9,angle(self.dataList{n}),self.plotSpecList{n});
                elseif strcmp(self.dataFormat,'degrees')
                    plot(self.frequencyList{n}/1e9,angle(self.dataList{n})*180/pi,self.plotSpecList{n});
                elseif strcmp(self.dataFormat,'magnitude')
                    plot(self.frequencyList{n}/1e9,abs(self.dataList{n}),self.plotSpecList{n});
                end
            end
            hold off
            xlabel('GHz')
            ylabel(self.dataFormat)
            title(self.title)
            legend(self.labelList,'Location','Best','Interpreter','None')
            grid on
        end
    end
end


% 	def addItem(self,data,arg2='S',arg3=(1,1),lineStyle=''):
% 		# (dataObj, dataItem, index)
% 		# (data, frequency, label)
% 		# interpret input arguments
% 		if isinstance(data,FrequencyDomainData):
% 			if isinstance(data,SParameters):
% 				self._dataList.append(eval("data."+arg2)[arg3[0]-1,arg3[1]-1,:])
% 				self._frequencyList.append(data.frequency)
% 				self.labelList.append(data.label + " " + arg2 + str(arg3))
% 				self.plotSpecList.append(lineStyle)
% 			elif isinstance(data,SpecLine):
% 				self._dataList.append(data.specLine)
% 				self._frequencyList.append(data.frequency)
% 				self.labelList.append(data.standard + " " + data.specItem)
% 				if lineStyle == '':
% 					lineStyle = '--'
% 				self.plotSpecList.append(lineStyle)
% 		elif type(data) == np.ndarray:
% 			self._dataList.append(data)
% 			self._frequencyList.append(arg2)
% 			self.labelList.append(arg3)
% 			self.plotSpecList.append(lineStyle)
% 
% 	def generatePlot(self):
% 		plt.figure()
% 		plt.hold(True)
% 		try:
% 			xlimits = self._xlimits
% 			if xlimits[0] == None:
% 				xlow = 0.
% 			else:
% 				xlow = xlimits[0]*1e9
% 			if xlimits[1] == None:
% 				xhigh = 1000e9
% 			else:
% 				xhigh = self._xlimits[1]*1e9
% 		except AttributeError:
% 			xlow = 0.
% 			xhigh = 1000e9
% 		for n in range(len(self._dataList)):
% 			if self.dataFormat == 'dB':
% 				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
% 					dB(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
% 					self.plotSpecList[n])
% 			elif self.dataFormat == 'radians':
% 				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
% 					np.angle(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
% 					self.plotSpecList[n])
% 			elif self.dataFormat == 'degrees':
% 				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
% 					180/np.pi*np.angle(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
% 					self.plotSpecList[n])
% 			elif self.dataFormat == 'magnitude':
% 				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
% 					np.abs(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
% 					self.plotSpecList[n])
% 		plt.hold(False)
% 		try:
% 			plt.xlim(self._xlimits)
% 		except (ValueError, AttributeError): 
% 			pass
% 		try:
% 			plt.ylim(self._ylimits)
% 		except (ValueError, AttributeError): 
% 			pass
% 		plt.xlabel('GHz')
% 		plt.ylabel(self.dataFormat)
% 		plt.title(self.title)
% 		plt.legend(self.labelList,fontsize='small')
% 		plt.grid(True)
