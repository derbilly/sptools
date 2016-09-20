# import some module functions
# numpy, matplotlib and scipy need to be installed
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import copy
import re
from scipy.optimize import curve_fit

def dB(val):
	return 20*np.log10(np.abs(val))

class FrequencyDomainData:
	def __init__(self):
		pass
	
class SParameters(FrequencyDomainData):
	def __init__(self,data):
		# data: touchstone file
		# data: SParameter object
		# data: dict with 'S', 'frequency' (,'portZ', 'label')
		FrequencyDomainData.__init__(self)
		if isinstance(data,SParameters):
			self.__S = data.S
			self.__frequency = data.frequency
			self.__numPorts = data.numPorts
			self.__portZ = data.portZ
			self.__dataFile = data.dataFile
			self.__label = data.label
		else:
			try: # data is touchstone file
				self.__dataFile = data
				self.__getSnp()
				# generate default label from filename
				pattern0 = re.compile('^.*/')
				pattern1 = re.compile('[.]s[0-9]*p')
				self.__label  = pattern1.sub('',pattern0.sub('',self.dataFile))
			except: 
				try: # data is dict
					self.__S = data['S']
					self.__frequency = data['frequency']
					try:
						self.__numPorts = data['numPorts']
					except KeyError:
						self.__numPorts = self.S.shape[1]
					try:
						self.__portZ = data['portZ']
					except KeyError:
						self.__portZ = 50.
					try:
						self.__dataFile = data['dataFile']
					except KeyError:
						self.__dataFile = ''
					try:
						self.__label = data['label']
					except KeyError:
						# generate default label from filename
						pattern0 = re.compile('^.*/')
						pattern1 = re.compile('[.]s[0-9]*p')
						self.__label  = pattern1.sub('',pattern0.sub('',self.dataFile))
				except:
					print("No such file or data: {0}".format(data))
					raise 
			
	def __str__(self):
		return ("### SParameter object ###\n" +
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
		return self.__S.copy()
	
	def getS(self,in1=(1,1),in2=None,frequency=None):
		"""getS(index) for index=(n,m) indexed from 1,
		returns Snm as a vector over the frequency index
		Usage:
			getS(n,m[,frequency])  indices as separate args with optional frequency
			getS((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.S[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	@property
	def frequency(self):
		return self.__frequency.copy()
		
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
	
	def plotAll(self,matrix='S'):
		plt.figure(figsize=(12,10))
		N = getattr(self,matrix).shape[0]
		for n in range(N):
			for m in range(N):
				plt.subplot(N,N,(n)*N+m+1)
				plt.plot(self.frequency*1e-9,dB(getattr(self,matrix)[n,m,:]))
				plt.xlabel('GHz')
				plt.ylabel('dB')
				plt.title(matrix+'('+str(n+1)+','+str(m+1)+')')
				plt.grid(True)
		plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)
		
	def reorderPorts(self,portOrder):
		"""portOrder is iterable of port indices indexed from 1"""
		####### check input
		portOrder = np.array(portOrder)-1
		self.__S = self.S[portOrder,:,:][:,portOrder,:]
		# set numPorts
		self.__numPorts = len(portOrder)
	
	def resampleFrequency(self,newFrequency):
		newS = np.zeros((self.numPorts,self.numPorts,len(newFrequency)),dtype=complex)
		for n in range(self.numPorts):
			for m in range(self.numPorts):
				newS[n,m,:] = fInterpolate(self.S[n,m,:],self.frequency,newFrequency)
		self.__S = newS
		self.__frequency = newFrequency
	
	def __getSnp(self):
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

		print("Importing "+self.dataFile)
		# open file for read
		fileID = open(self.dataFile,'r')
		# read file
		# read header: # frequencyUnit networkFormat complexFormat R portZ (complex impedance terms not allowed)
		nextLine = getNextLine(fileID) # header line ### error check
		header = nextLine.split()
		frequencyUnit = header[1].upper()
		networkFormat = header[2].upper()
		complexFormat = header[3].upper()
		self.__portZ = float(header[5])
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
		self.__numPorts = int(np.sqrt(numComplexValuesPerFrequencyPoint)) ### error check this
		# read and create data matrices
		self.__frequency = np.zeros(numFrequencyPoints)
		networkData1 = np.zeros((self.numPorts,self.numPorts,numFrequencyPoints))
		networkData2 = np.zeros((self.numPorts,self.numPorts,numFrequencyPoints))
		for nf in range(numFrequencyPoints):
			self.__frequency[nf] = float(lines[nf*numLinesPerFrequencyPoint].split()[0])
			dataItems = lines[nf*numLinesPerFrequencyPoint].split()[1:]
			for nlpfp in range(1,numLinesPerFrequencyPoint):
				dataItems.extend(lines[nf*numLinesPerFrequencyPoint+nlpfp].split())
			for n in range(self.numPorts): # row
				for m in range(self.numPorts): # column
					networkData1[n,m,nf] = float(dataItems[(n*self.numPorts+m)*2])
					networkData2[n,m,nf] = float(dataItems[(n*self.numPorts+m)*2+1])
		# scale frequency values to Hz
		if frequencyUnit == 'HZ':
			frequencyFactor = 1e0
		elif frequencyUnit == 'KHZ':
			frequencyFactor = 1e3
		elif frequencyUnit == 'MHZ':
			frequencyFactor = 1e6
		elif frequencyUnit == 'GHZ':
			frequencyFactor = 1e9
		self.__frequency *= frequencyFactor
		# generate complex valued data from input data format
		self.__S = np.zeros((self.numPorts,self.numPorts,numFrequencyPoints),dtype=complex)
		if complexFormat == 'RI':
			self.__S = networkData1 + networkData2*1j
		elif complexFormat == 'MA':
			self.__S = networkData1 * np.exp(1j*2*np.pi/360.0*networkData2)
		elif complexFormat == 'DB':
			self.__S = 10**(networkData1/20.0) * np.exp(1j*2*np.pi/360.0*networkData2)
	
	def export(self, dataFile='export', dataFormat='RI', numDigits=12, format='touchstone', frequencyFormat='GHZ'):
		# remove snp file extension and add a proper one
		regexp0 = re.compile(r"\.s[0-9]*p$")
		dataFile = regexp0.sub("",dataFile)
		dataFile = "{0}.s{1}p".format(dataFile,self.numPorts)
		
		print("Exporting " + dataFile)
		
		dataFormat = dataFormat.upper()
		if dataFormat == 'RI':
			data1 = np.real(self.S)
			data2 = np.imag(self.S)
		elif dataFormat == 'MA':
			data1 = np.abs(self.S)
			data2 = np.angle(self.S)
		elif dataFormat == 'DB':
			data1 = dB(self.S)
			data2 = np.angle(self.S)
		else:
			dataFormat = 'RI'
			data1 = np.real(self.S)
			data2 = np.imag(self.S)
			
		freq = self.frequency*1e-9
		fid = open(dataFile,'w')
		regexp1 = re.compile(r"\n") # add ! comment character to lines
		fid.write("! Touchstone data exported from sptools\n! " +
			regexp1.sub("\n! ",str(self)) )
		fid.write("\n# GHZ S {0} R {1}".format(dataFormat,self.portZ))
		for frequencyIndex in range(len(freq)):
			fid.write("\n{0:0.5f}\t".format(freq[frequencyIndex]))
			for n in range(self.numPorts):
				for m in range(self.numPorts):
					if (n*self.numPorts + m) != 0 and (n*self.numPorts + m) % 4 == 0:
						fid.write("\n\t\t")
					fid.write("{0:0.14g} {1:0.14g} ".format(data1[n,m,frequencyIndex],data2[n,m,frequencyIndex]))
		fid.write("\n")
		fid.close()

class MixedModeSParameters(SParameters):		
	def __init__(self,data):
		SParameters.__init__(self,data)
		self.__genSMM()
	
	def __str__(self):
		return ("### MixedModeSParameter object ###\n" +
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
		return self.__SMM.copy()
	
	def getSMM(self,in1=(1,1),in2=None,frequency=None):
		"""getSMM(index) for index=(n,m) indexed from 1,
		returns SMMnm as a vector over the frequency index
		Usage:
			getSMM(n,m[,frequency])  indices as separate args with optional frequency
			getSMM((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.SMM[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	@property
	def SDD(self):
		return self.SMM[0::2,0::2]
	
	def getSDD(self,in1=(1,1),in2=None,frequency=None):
		"""getSDD(index) for index=(n,m) indexed from 1,
		returns SDDnm as a vector over the frequency index
		Usage:
			getSDD(n,m[,frequency])  indices as separate args with optional frequency
			getSDD((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.SDD[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	@property
	def SDC(self):
		return self.SMM[0::2,1::2]
	
	def getSDC(self,in1=(1,1),in2=None,frequency=None):
		"""getSDC(index) for index=(n,m) indexed from 1,
		returns SDCnm as a vector over the frequency index
		Usage:
			getSDC(n,m[,frequency])  indices as separate args with optional frequency
			getSDC((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.SDC[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	@property
	def SCD(self):
		return self.SMM[1::2,0::2]
	
	def getSCD(self,in1=(1,1),in2=None,frequency=None):
		"""getSCD(index) for index=(n,m) indexed from 1,
		returns SCDnm as a vector over the frequency index
		Usage:
			getSCD(n,m[,frequency])  indices as separate args with optional frequency
			getSCD((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.SCD[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	@property
	def SCC(self):
		return self.SMM[1::2,1::2]
	
	def getSCC(self,in1=(1,1),in2=None,frequency=None):
		"""getSCC(index) for index=(n,m) indexed from 1,
		returns SCCnm as a vector over the frequency index
		Usage:
			getSCC(n,m[,frequency])  indices as separate args with optional frequency
			getSCC((n,m)[,frequency]) indices as a tuple with optional frequency
			frequency can be specified as a single value or list of values
		"""
		try:
			n = in1[0]-1
			m = in1[1]-1
			if frequency == None:
				frequency = in2
		except TypeError:
			n = in1 - 1
			m = in2 - 1

		data = self.SCC[n,m,:]
		freq = self.frequency

		if np.max(frequency)==None:
			return data
		else:
			freqreq = np.array(frequency) # requested frequencies
			if np.max(freqreq) <= np.max(freq) and np.min(freqreq) >= np.min(freq):
				if freqreq in freq: # 
					return data[freq==freqreq]
				else: # interpolate the values
					return fInterpolate(data,freq,freqreq)
			else:
				print("Frequency out of range")

	# other functions
	def copy(self):
		return copy.deepcopy(self)

	def __genSMM(self):
		if np.size(self.S,0)%2 != 0:
			raise MixedModeOddNumberOfPorts
		numPorts = np.size(self.S,0)//2 # mixed mode ports
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
		self.__SMM = np.zeros_like(self.S)
		for n in range(np.size(self.S,2)):
			self.__SMM[:,:,n] = np.dot( M, np.dot(self.S[:,:,n],Mi) ) # M*S*Mi
	
	def reorderPorts(self,portOrder):
		SParameters.reorderPorts(self,portOrder)
		self.__genSMM()
	
	def resampleFrequency(self,newFrequency):
		SParameters.resampleFrequency(self,newFrequency)
		self.__genSMM()

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
	
	# 10GBASE-KR
	# IEEE 802.3 Annex69B
	
	def gen_10GKR_Amax(self):
		b1 = 2.0e-5
		b2 = 1.1e-10
		b3 = 3.2e-20
		b4 = -1.2e-30
		f1 = 1.0e9
		f2 = 6.0e9
		self.frequency = np.linspace(f1, f2, int((f2-f1)/10e6)+1)
		dBspecLine = -20*np.log10(np.exp(1))*(b1*np.sqrt(self.frequency) + b2*self.frequency +
			b3*self.frequency**2 + b4*self.frequency**3)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_10GKR_ILmax(self):
		b1 = 2.0e-5
		b2 = 1.1e-10
		b3 = 3.2e-20
		b4 = -1.2e-30
		f2 = 6.0e9
		fmin = 0.05e9
		fmax = 15.0e9
		self.frequency = np.linspace(fmin, fmax, int((fmax-fmin)/10e6)+1)
		dBspecLine = -20*np.log10(np.exp(1))*(b1*np.sqrt(self.frequency) + b2*self.frequency +
			b3*self.frequency**2 + b4*self.frequency**3) - 0.8 - 2.0e-10*self.frequency
		dBspecLine[self.frequency>f2] = dBspecLine[self.frequency>f2] - 1e-8*(self.frequency[self.frequency>f2]-f2)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_10GKR_ILDmin(self):
		f1 = 1.0e9
		f2 = 6.0e9
		self.frequency = np.linspace(f1, f2, int((f2-f1)/10e6)+1)
		dBspecLine = -(1.0 + 0.5e-9*self.frequency)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'

	def gen_10GKR_ILDmax(self):
		f1 = 1.0e9
		f2 = 6.0e9
		self.frequency = np.linspace(f1, f2, int((f2-f1)/10e6)+1)
		dBspecLine = (1.0 + 0.5e-9*self.frequency)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'

	def gen_10GKR_RLmin(self):
		f2 = 10312.5e6
		self.frequency = np.linspace(50e6, f2, int((f2-50e6)/12.5e6)+1)
		dBspecLine = -(12.0 - 6.75*np.log10(self.frequency/275e6))
		dBspecLine[self.frequency<275.0e6] = -12
		dBspecLine[self.frequency>3000.0e6] = -5
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'

	def gen_10GKR_ICRmin(self):
		fa = 0.1e9
		fb = 5.15625e9
		self.frequency = np.linspace(fa, fb, 501)
		dBspecLine = 23.3 -  18.7*np.log10(self.frequency/5.0e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'


	generatorMethod['10GBASE-KR'] = {
	'Amax':gen_10GKR_Amax, # 69B.4.2 Fitted attenuation
	'ILmax':gen_10GKR_ILmax, # 69B.4.3 Insertion loss
	'ILDmin':gen_10GKR_ILDmin, # 69B.4.4 Insertion loss deviation
	'ILDmax':gen_10GKR_ILDmax, # 69B.4.4 Insertion loss deviation
	'RLmin':gen_10GKR_RLmin, # 69B.4.5 Return loss
	'ICRmin':gen_10GKR_ICRmin} # 69B.4.6 Crosstalk

	
	
	# chip-to-module CAUI4
	# IEEE 802.3bm
	def gen_ctmCAUI4_IL(self):
		self.frequency = np.linspace(10e6,18.75e9,int(18.75e9/10e6))
		dBspecLine = -(1.076*(-18 + 2*self.frequency/1e9))
		dBspecLine[self.frequency<14e9] = -(1.076*(0.075 +
			0.537*np.sqrt(self.frequency[self.frequency<14e9]/1e9) + 
			0.566*self.frequency[self.frequency<14e9]/1e9))
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_ctmCAUI4_RLd(self):
		self.frequency = np.linspace(10e6,19e9,int(19e9/10e6))
		dBspecLine = -(4.75-7.4*np.log10(self.frequency/14e9))
		dBspecLine[self.frequency<8e9] = -(9.5 - 
			0.37*self.frequency[self.frequency<8e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'
	
	def gen_ctmCAUI4_RLdc(self):
		self.frequency = np.linspace(10e6,19e9,int(19e9/10e6))
		dBspecLine = -(15-6*self.frequency/25.78e9)
		dBspecLine[self.frequency<12.89e9] = -(22 - 
			20*self.frequency[self.frequency<12.89e9]/1e9/25.78)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'
	
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
		self.limitType = 'min'
	
	def gen_100GCR4_IL_catf(self): # 92.11.2 Cable assembly test fixture
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(-0.00125 + 0.12*np.sqrt(self.frequency/1e9) + 0.0575*self.frequency/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_100GCR4_IL_MTFmax(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(-4.25 + 0.66*self.frequency/1e9)
		dBspecLine[self.frequency<14e9] = -(0.12 +
			0.475*np.sqrt(self.frequency[self.frequency<14e9]/1e9) + 
			0.221*self.frequency[self.frequency<14e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_100GCR4_IL_MTFmin(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(0.0656*np.sqrt(self.frequency/1e9) + 0.164*self.frequency/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'
	
	def gen_100GCR4_RLd_MTF(self): # 92.11.3 Mated test fixtures
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(18 - 0.5*self.frequency/1e9)
		dBspecLine[self.frequency<4e9] = -(20 -
			1*self.frequency[self.frequency<4e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'
	
	def gen_100GCR4_ILdc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(30 - 29/22*self.frequency/1e9)
		dBspecLine[self.frequency>=16.5e9] = -8.25
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'min'
	
	def gen_100GCR4_RLc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(12 - 9*self.frequency/1e9)
		dBspecLine[self.frequency>=1e9] = -3
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'

	def gen_100GCR4_RLdc_MTF(self): # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		self.frequency = np.linspace(10e6,25e9,int(25e9/10e6))
		dBspecLine = -(18 - 6/25.78*self.frequency/1e9)
		dBspecLine[self.frequency<12.89e9] = -(30 -
			30/25.78*self.frequency[self.frequency<12.89e9]/1e9)
		self.specLine = 10**(dBspecLine/20)
		self.limitType = 'max'
	
	generatorMethod['100GBASE-CR4'] = {
		'IL_tfref':gen_100GCR4_IL_tfref, # 92.11.1.2 Test fixture insertion loss
		'IL_catf':gen_100GCR4_IL_catf, # 92.11.2 Cable assembly test fixture
		'IL_MTFmax':gen_100GCR4_IL_MTFmax, # 92.11.3 Mated test fixtures
		'IL_MTFmin':gen_100GCR4_IL_MTFmin, # 92.11.3 Mated test fixtures
		'RLd_MTF':gen_100GCR4_RLd_MTF, # 92.11.3 Mated test fixtures
		'ILdc_MTF':gen_100GCR4_ILdc_MTF, # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		'RLc_MTF':gen_100GCR4_RLc_MTF, # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss
		'RLdc_MTF':gen_100GCR4_RLdc_MTF} # 92.11.3.3 Mated test fixtures common-mode conversion insertion loss

class DataPlot:
	def __init__(self):
		pass
		
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
		DataPlot.__init__(self)
		self.dataFormat = dataFormat
		self.title = title
		self._dataList = []
		self._frequencyList = []
		self.labelList = []
		self.plotSpecList = []
		
	def addItem(self,data,arg2='S',arg3=(1,1),lineStyle=''):
		# (dataObj, dataItem, index)
		# (data, frequency, label)
		# interpret input arguments
		if isinstance(data,FrequencyDomainData):
			if isinstance(data,SParameters):
				self._dataList.append(eval("data."+arg2)[arg3[0]-1,arg3[1]-1,:])
				self._frequencyList.append(data.frequency)
				self.labelList.append(data.label + " " + arg2 + str(arg3))
				self.plotSpecList.append(lineStyle)
			elif isinstance(data,SpecLine):
				self._dataList.append(data.specLine)
				self._frequencyList.append(data.frequency)
				self.labelList.append(data.standard + " " + data.specItem)
				if lineStyle == '':
					lineStyle = '--'
				self.plotSpecList.append(lineStyle)
		elif type(data) == np.ndarray:
			self._dataList.append(data)
			self._frequencyList.append(arg2)
			self.labelList.append(arg3)
			self.plotSpecList.append(lineStyle)

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
		
def calculateLoss(data):
	(numPorts,_,numFrequencyPoints)=data.shape
	loss = np.ones((numPorts,numFrequencyPoints))
	power = np.zeros((numPorts,numFrequencyPoints))
	for m in range(numPorts):
		for n in range(numPorts):
			loss[m,:] -= np.abs(data[n,m,:])**2
			power[m,:] += np.abs(data[n,m,:])**2
	return np.sqrt(power)

def fInterpolate(data,f,fnew):
	magnitude = np.abs(data)
	phase = np.unwrap(np.angle(data))
	newMagnitude = interpolate.pchip_interpolate(f,magnitude,fnew)
	newPhase = interpolate.pchip_interpolate(f,phase,fnew)
	newData = newMagnitude*np.exp(1j*newPhase)
	return newData

class BackplaneEthernetChannel:
	# Annex 69B
	
	# constants
	f_min = 0.05e9
	f_max = 15.00e9
	def __init__(self,standard='10GBASE-KR'):
		self.__thru = {}
		self.__returnLoss = []
		self.__NEXT = []
		self.__FEXT = []
		if standard=='10GBASE-KR':
			self.f_1 = 1.0e9
			self.f_2 = 6.0e9
			self.f_a = 0.1e9
			self.f_b = 5.15625e9
	@property
	def thru(self):
		return self.__thru.copy()
	
	def addThru(self,data,frequency,index=(2,1)):	# with 2 port S params or thru, add RL if available
		if len(data.shape) == 1:
			self.__thru = {'data':data,'frequency':frequency}
		elif len(data.shape) == 3:
			self.__thru = {'data':data[1,0,:],'frequency':frequency}
			self.__returnLoss.append({'data':data[0,0,:],'frequency':frequency})
			self.__returnLoss.append({'data':data[1,1,:],'frequency':frequency})
	
	def addRL(self,data,frequency):	# add RL
		self.__returnLoss.append({'data':data,'frequency':frequency})
	
	def addNEXT(self,data,frequency):
		if len(data.shape) == 1:
			self.__NEXT.append({'data':data,'frequency':frequency})
		elif len(data.shape) == 3:
			self.__NEXT.append({'data':data[1,0,:],'frequency':frequency})
	
	def addFEXT(self,data,frequency):
		if len(data.shape) == 1:
			self.__FEXT.append({'data':data,'frequency':frequency})
		elif len(data.shape) == 3:
			self.__FEXT.append({'data':data[1,0,:],'frequency':frequency})
	
	def __calcILD(self):
		fi = np.linspace(self.f_1, self.f_2, 501)
		IL = -dB(fInterpolate(self.__thru['data'], self.__thru['frequency'], fi))
		f_avg = np.sum(fi)/len(fi)
		IL_avg = np.sum(IL)/len(IL)
		m_A = sum((fi-f_avg)*(IL-IL_avg)) / sum((fi-f_avg)**2)
		b_A = IL_avg - m_A*f_avg
		AdB = -(m_A*self.__thru['frequency'] + b_A)
		self.__thru['A'] = 10**(AdB/20)
		self.__thru['A'][self.__thru['frequency'] < self.f_1] = np.nan
		self.__thru['A'][self.__thru['frequency'] > self.f_2] = np.nan
		self.__thru['ILD'] = 10**(-(dB(self.thru['data']) - dB(self.thru['A']))/20)
		
	def __calcXT(self):
		# need to interpolate the data if necessary
		# all PS*XT expressed as mag: use dB to represent
		# all quantities
		xtsum = np.zeros(len(self.thru['frequency']))
		for xt in self.__NEXT:
			xtsum += 10**(dB(xt['data'])/10)
		self.PSNEXT = 10**(10*np.log10(xtsum)/20)
		xtsum = np.zeros(len(self.thru['frequency']))
		for xt in self.__FEXT:
			xtsum += 10**(dB(xt['data'])/10)
		self.PSFEXT = 10**(10*np.log10(xtsum)/20)
		self.PSXT = 10**(10*np.log10(10**(dB(self.PSNEXT)/10) + 10**(dB(self.PSFEXT)/10))/20)
		
	def __calcICR(self):
		ICRdB = dB(self.thru['data']) - dB(self.PSXT)
		fi = np.linspace(self.f_a, self.f_b, 501)
		f_avg = np.sum(fi)/len(fi)
		#ICRavg = 1/
	# fitted Attenuation
	# IL
	# RL
	# ILD
	# PSNEXT
	# PSFEXT
	# PSXT
	# ICR
	
	def evaluate(self):
		self.__calcILD()
		self.__calcXT()
		self.__calcICR()
		
def limitMargin(data,limit,dataFormat='dB'):
	# limitMargin(data,limit,dataFormat='dB')
	# data: (data,freq) e.g. (x.getSDD(2,1),x.frequency)
	# limit: SpecLine
	# returns tuple of (margin, frequency)
	pass