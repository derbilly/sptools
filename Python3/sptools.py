# import some module functions
# numpy, matplotlib and scipy need to be installed
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import copy
import re

def dB(val):
	return 20*np.log10(np.abs(val))

def getSnp(dataFile):
	# local functions
	def getNextLine(fileID):
		while True:
			# read line with \n
			line = fileID.readline()
			# check if it is EOF, ''
			if line=='':
				return -1
			elif not(isComment(line)):
				# return line without newline
				return line[0:-1]

	def isComment(line):
		# remove newline and spaces
		line = line[0:-1].replace(' ','')
		if not(bool(line)):
			return True
		elif line[0]=='!':
			return True
		else:
			return False


	# open file for read
	fileID = open(dataFile,'r')
	# read file
	# read header: # frequencyUnit networkFormat complexFormat R portZ (complex impedance terms not allowed)
	nextLine = getNextLine(fileID) # header line ### error check
	header = nextLine.split()
	frequencyUnit = header[1].upper()
	networkFormat = header[2].upper()
	complexFormat = header[3].upper()
	portZ = float(header[5])
	#### check for errors in header line
	# read data
	# read all remaining lines
	lines = []
	nextLine = getNextLine(fileID)
	while nextLine != -1:
		lines.append(nextLine)
		nextLine = getNextLine(fileID)
	# close file
	fileID.close()
	# find how many lines per frequency point
	numLinesPerFrequencyPoint = 1
	numComplexValuesPerFrequencyPoint = 0
	if len(lines[0].split())%2 == 0:
		raise touchstoneFormatError # should be odd with freq value included with complex data
	numComplexValuesPerFrequencyPoint += (len(lines[0].split()) - 1)/2
	while len(lines[numLinesPerFrequencyPoint].split())%2 == 0:
		numLinesPerFrequencyPoint += 1
		numComplexValuesPerFrequencyPoint += len(lines[numLinesPerFrequencyPoint].split())/2
	# check that there are an expected number of total lines, i.e. numLinesPerFrequencyPoint * number of frequency points
	if len(lines)%numLinesPerFrequencyPoint != 0:
		raise touchstoneFormatError # some extra junk at end
	numFrequencyPoints = len(lines)//numLinesPerFrequencyPoint
	numPorts = int(np.sqrt(numComplexValuesPerFrequencyPoint)) ### error check this
	# read and create data matrices
	frequency = np.zeros(numFrequencyPoints)
	networkData1 = np.zeros((numPorts,numPorts,numFrequencyPoints))
	networkData2 = np.zeros((numPorts,numPorts,numFrequencyPoints))
	for nf in range(numFrequencyPoints):
		frequency[nf] = float(lines[nf*numLinesPerFrequencyPoint].split()[0])
		dataItems = lines[nf*numLinesPerFrequencyPoint].split()[1:]
		for nlpfp in range(1,numLinesPerFrequencyPoint):
			dataItems.extend(lines[nf*numLinesPerFrequencyPoint+nlpfp].split())
		for n in range(numPorts): # row
			for m in range(numPorts): # column
				networkData1[n,m,nf] = float(dataItems[(n*numPorts+m)*2])
				networkData2[n,m,nf] = float(dataItems[(n*numPorts+m)*2+1])
	# scale frequency values to Hz
	if frequencyUnit == 'HZ':
		frequencyFactor = 1e0
	elif frequencyUnit == 'KHZ':
		frequencyFactor = 1e3
	elif frequencyUnit == 'MHZ':
		frequencyFactor = 1e6
	elif frequencyUnit == 'GHZ':
		frequencyFactor = 1e9
	frequency *= frequencyFactor
	# generate complex valued data from input data format
	networkData = np.zeros((numPorts,numPorts,numFrequencyPoints),dtype=complex)
	if complexFormat == 'RI':
		networkData = networkData1 + networkData2*1j
	elif complexFormat == 'MA':
		networkData = networkData1 * np.exp(1j*2*np.pi/360.0*networkData2)
	elif complexFormat == 'DB':
		networkData = 10**(networkData1/20.0) * np.exp(1j*2*np.pi/360.0*networkData2)
	# return data
	return (networkData, frequency, numPorts, portZ)

def genSMM(S):
	if np.size(S,0)%2 != 0:
		raise MixedModeOddNumberOfPorts
	numPorts = np.size(S,0)//2 # mixed mode ports
	# construct transformation matrix
	if numPorts == 1: # add case for single mm port
		M = np.array([[1,-1],[1,1]])/np.sqrt(2)
	else:
		M0 = np.array([[1,0,-1],[1,0,1]])/np.sqrt(2)
		M = np.zeros((2*numPorts,2*numPorts))
		for n in range(0,numPorts,2):
			M[  2*n:2*n+2 ,   2*n:2*n+3] = M0
			M[2*n+2:2*n+4 , 2*n+1:2*n+4] = M0
	Mi = np.linalg.inv(M)
	SMM = np.zeros_like(S)
	for n in range(np.size(S,2)):
		SMM[:,:,n] = np.dot( M, np.dot(S[:,:,n],Mi) ) # M*S*Mi
	return SMM
	
class FrequencyDomainData:
	def __init__(self):
		pass
	
class SParameters(FrequencyDomainData):
	def __init__(self,dataFile):
		FrequencyDomainData.__init__(self)
		self.__dataFile = dataFile
		(self.__S, self.__frequency, self.__numPorts, self.__portZ)  = getSnp(self.__dataFile)
		pattern0 = re.compile('^.*/')
		pattern1 = re.compile('[.]s[0-9]*p')
		self.__label  = pattern1.sub('',pattern0.sub('',dataFile))
			
	# @classmethod
	# def fromTouchstoneFile(cls, dataFile):
		# self.__dataFile = dataFile
		# (self.__S, self.__frequency, self.__numPorts, self.__portZ)  = getSnp(self.__dataFile)
		# pattern0 = re.compile('^.*/')
		# pattern1 = re.compile('[.]s[0-9]*p')
		# self.__label  = pattern1.sub('',pattern0.sub('',dataFile))
	
	
	def __str__(self):
		return ("SParameter object\n" +
				"Label:           " + self.label + "\n" +
				"Datafile:        " + self.dataFile + "\n" +
				"Number of ports: " + str(self.numPorts) + "\n" +
				"Port impedance:  " + str(self.portZ) + " Ohms\n" +
				"Frequency start: " + str(self.frequency[0]*1e-9) + " GHz\n" +
				"Frequency stop:  " + str(self.frequency[-1]*1e-9) + " GHz\n" +
				"Frequency steps: " + str(len(self.frequency)))
	
	# setters and getters
	@property
	def label(self):
		return self.__label
	
	@label.setter
	def label(self,label):
		self.__label = label
		
	@property
	def S(self):
		"""S-Parameters:
		S(row,col,freq)	- ndarray
			row:	destination port, index from 0
			col:	source port, index from 0
			freq:	frequency index, index from 0
		"""
		return self.__S
	
	def getS(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.S[n,m,:]

	@property
	def frequency(self):
		return self.__frequency
		
	@property
	def numPorts(self):
		return self.__numPorts
	
	@property
	def portZ(self):
		return self.__portZ
	
	@property
	def dataFile(self):
		return self.__dataFile
	
	# more functions
	def copy(self):
		return copy.deepcopy(self)
		
	def reorderPorts(self,portOrder):
		"""portOrder is iterable of port indices indexed from 1"""
		####### check input
		portOrder = np.array(portOrder)-1
		self.__S = self.S[portOrder,:,:][:,portOrder,:]
		# set numPorts
		self.__numPorts = len(portOrder)
	
	def resampleFrequency(self,newFrequency):
		newS = np.zeros((self.numPorts,self.numPorts,len(newFrequency)))
		for n in range(self.numPorts):
			for m in range(self.numPorts):
				magnitude = np.abs(self.S[n,m,:])
				phase = np.unwrap(np.angle(self.S[n,m,:]))
				newMagnitude = interpolate.pchip_interpolate(self.frequency,magnitude,newFrequency)
				newPhase = interpolate.pchip_interpolate(self.frequency,phase,newFrequency)
				newS[n,m,:] = newMagnitude*np.exp(1j*newPhase)
		self.__S = newS
		self.__frequency = newFrequency
		
	def exportTouchstoneFile(self, dataFile='export', dataFormat='RI', numSignificantDigits=8):
		pass

class MixedModeSParameters(SParameters):		
	def __init__(self,touchstoneFile):
		SParameters.__init__(self,touchstoneFile)
		self.__SMM = genSMM(self.S)
	
	def __str__(self):
		return ("MixedModeSParameter object\n" +
				"Label:           " + self.label + "\n" +
				"Datafile:        " + self.dataFile + "\n" +
				"Number of ports: " + str(self.numPorts) + "\n" +
				"Port impedance:  " + str(self.portZ) + " Ohms\n" +
				"Frequency start: " + str(self.frequency[0]*1e-9) + " GHz\n" +
				"Frequency stop:  " + str(self.frequency[-1]*1e-9) + " GHz\n" +
				"Frequency steps: " + str(len(self.frequency)))
	
	# setters and getters
	@property
	def SMM(self):
		return self.__SMM
	
	def getSMM(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.SMM[n,m,:]

	@property
	def SDD(self):
		return self.SMM[0::2,0::2]
	
	def getSDD(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.SDD[n,m,:]

	@property
	def SDC(self):
		return self.SMM[0::2,1::2]
	
	def getSDC(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.SDC[n,m,:]

	@property
	def SCD(self):
		return self.SMM[1::2,0::2]
	
	def getSCD(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.SCD[n,m,:]

	@property
	def SCC(self):
		return self.SMM[1::2,1::2]
	
	def getSCC(self,n0=1,m0=1):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		"""
		if type(n0) == int:
			n=n0-1
			m=m0-1
		elif len(n0)==2:
			n = n0[0]-1
			m = n0[1]-1
		return self.SCC[n,m,:]

	# other functions
	def copy(self):
		return copy.deepcopy(self)

	def reorderPorts(self,portOrder):
		SParameters.reorderPorts(self,portOrder)
		self.__SMM = genSMM(self.S)
	
	def resampleFrequency(self,newFrequency):
		SParameters.resampleFrequency(self,newFrequency)
		self.__SMM = genSMM(self.S)

class SpecLine(FrequencyDomainData):
	def __init__(self,standard,specItem):
		FrequencyDomainData.__init__(self)
		self.specItem = specItem
		self.standard = standard
		try:
			self.generatorMethod[standard][specItem](self)
		except KeyError:
			try:
				subdict = self.generatorMethod[standard]
				print("Spec item: '" + self.specItem + "' not implemented.")
				print("Valid items:")
				print("\t"+'\n\t'.join(subdict.keys()))
			except KeyError:
				print("Standard: '" + self.standard + "' not implemented.")
				print("Valid standards:")
				print ("\t"+'\n\t'.join(self.generatorMethod.keys()))
				
		
	# SpecLine generator methods dictionary
	generatorMethod = dict()

	# chip-to-module CAUI4
	# IEEE 802.3bm
	def gen_ctmCAUI4_IL(self):
		self.frequency = np.linspace(10e6,18.75e9,int(18.75e9/10e6))
		dBspecLine = -(1.076*(-18 + 2*self.frequency/1e9))
		dBspecLine[self.frequency<14e9] = -(1.076*(0.075 +
			0.537*np.sqrt(self.frequency[self.frequency<14e9]/1e9) + 
			0.566*self.frequency[self.frequency<14e9]/1e9))
		self.specLine = 10**(dBspecLine/20)
	
	def gen_ctmCAUI4_RLd(self):
		self.frequency = np.linspace(10e6,19e9,int(19e9/10e6))
		dBspecLine = -(4.75-7.4*np.log10(self.frequency/14e9))
		dBspecLine[self.frequency<8e9] = -(9.5 - 
			0.37*self.frequency[self.frequency<8e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_ctmCAUI4_RLdc(self):
		self.frequency = np.linspace(10e6,19e9,int(19e9/10e6))
		dBspecLine = -(15-6*self.frequency/25.78e9)
		dBspecLine[self.frequency<12.89e9] = -(22 - 
			20*self.frequency[self.frequency<12.89e9]/1e9/25.78)
		self.specLine = 10**(dBspecLine/20)
	
	generatorMethod['chip-to-module CAUI4'] = {
		'IL':gen_ctmCAUI4_IL,
		'RLd':gen_ctmCAUI4_RLd,
		'RLdc':gen_ctmCAUI4_RLdc}

	# 100GBASE-CR4
	# IEEE 802.3bj 92.11
	# Also used for chip-to-module CAUI4
	def gen_100GCR4_IL_tfref(self): # 92.11.1.2 Test fixture insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(-0.00144 + 0.13824*np.sqrt(self.frequency/1e9) + 0.06624*self.frequency/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_IL_catf(self): # 92.11.2 Cable assembly test fixture
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(-0.00125 + 0.12*np.sqrt(self.frequency/1e9) + 0.0575*self.frequency/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_IL_MTFmax(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(-4.25 + 0.66*self.frequency/1e9)
		dBspecLine[self.frequency<14e9] = -(0.12 +
			0.475*np.sqrt(self.frequency[self.frequency<14e9]/1e9) + 
			0.221*self.frequency[self.frequency<14e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_IL_MTFmin(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(0.0656*np.sqrt(self.frequency/1e9) + 0.164*self.frequency/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_RLd_MTF(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(18 - 0.5*self.frequency/1e9)
		dBspecLine[self.frequency<4e9] = -(20 -
			1*self.frequency[self.frequency<4e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_ILdc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(30 - 29/22*self.frequency/1e9)
		dBspecLine[self.frequency>=16.5e9] = -8.25
		self.specLine = 10**(dBspecLine/20)
	
	def gen_100GCR4_RLc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(12 - 9*self.frequency/1e9)
		dBspecLine[self.frequency>=1e9] = -3
		self.specLine = 10**(dBspecLine/20)

	def gen_100GCR4_RLdc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(18 - 6/25.78*self.frequency/1e9)
		dBspecLine[self.frequency<12.89e9] = -(30 -
			30/25.78*self.frequency[self.frequency<12.89e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
	
	generatorMethod['100GBASE-CR4'] = {
		'IL_tfref':gen_100GCR4_IL_tfref,
		'IL_catf':gen_100GCR4_IL_catf,
		'IL_MTFmax':gen_100GCR4_IL_MTFmax,
		'IL_MTFmin':gen_100GCR4_IL_MTFmin,
		'RLd_MTF':gen_100GCR4_RLd_MTF,
		'ILdc_MTF':gen_100GCR4_ILdc_MTF,
		'RLc_MTF':gen_100GCR4_RLc_MTF,
		'RLdc_MTF':gen_100GCR4_RLdc_MTF}

class DataPlot:
	def xlim(self,low=None,high=None):
		self._xlimits = (low,high)
	
	def ylim(self,low=None,high=None):
		self._ylimits = (low,high)
	
class FrequencyDomainPlot(DataPlot):
	# data format: dB, magnitude, phase, degrees
	# xlimits
	# ylimits reset for change type
	# title
	# list of data items, line formats and legend labels
	def __init__(self,dataFormat='dB',title=''):
		self.dataFormat = dataFormat
		self.title = title
		self._dataList = []
		self._frequencyList = []
		self.labelList = []
		self.plotSpecList = []
		
	def addItem(self,data,arg2='S',arg3=(1,1)):
		# (dataObj, dataItem, index)
		# (data, frequency, label)
		# interpret input arguments
		if isinstance(data,FrequencyDomainData):
			if isinstance(data,SParameters):
				self._dataList.append(eval("data."+arg2)[arg3[0]-1,arg3[1]-1,:])
				self._frequencyList.append(data.frequency)
				self.labelList.append(data.label + " " + arg2 + str(arg3))
				self.plotSpecList.append('')
			elif isinstance(data,SpecLine):
				self._dataList.append(data.specLine)
				self._frequencyList.append(data.frequency)
				self.labelList.append(data.standard + " " + data.specItem)
				self.plotSpecList.append('--')
		elif type(data) == np.ndarray:
			self._dataList.append(data)
			self._frequencyList.append(arg2)
			self.labelList.append(arg3)
			self.plotSpecList.append('')

	def generatePlot(self):
		plt.figure()
		plt.hold(True)
		try:
			xlimits = self._xlimits
			if xlimits[0] == None:
				xlow = 0.
			else:
				xlow = xlimits[0]*1e9
			if xlimits[1] == None:
				xhigh = 1000e9
			else:
				xhigh = self._xlimits[1]*1e9
		except AttributeError:
			xlow = 0.
			xhigh = 1000e9
		for n in range(len(self._dataList)):
			if self.dataFormat == 'dB':
				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
					dB(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
					self.plotSpecList[n])
			elif self.dataFormat == 'radians':
				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
					np.angle(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
					self.plotSpecList[n])
			elif self.dataFormat == 'degrees':
				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
					180/np.pi*np.angle(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
					self.plotSpecList[n])
			elif self.dataFormat == 'magnitude':
				plt.plot(self._frequencyList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]*1e-9,
					np.abs(self._dataList[n][np.logical_and(self._frequencyList[n]>=xlow, self._frequencyList[n]<=xhigh)]),
					self.plotSpecList[n])
		plt.hold(False)
		try:
			plt.xlim(self._xlimits)
		except (ValueError, AttributeError): 
			pass
		try:
			plt.ylim(self._ylimits)
		except (ValueError, AttributeError): 
			pass
		plt.xlabel('GHz')
		plt.ylabel(self.dataFormat)
		plt.title(self.title)
		plt.legend(self.labelList,fontsize='small')
		plt.grid(True)
		
