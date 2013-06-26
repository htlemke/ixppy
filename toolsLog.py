import utilities

g_loglevel=1

class Message(object):
	def __init__(self,msg,fromfunction=None):
		self.msg=msg
		self.time = utilities.now()
		self.fromfunction=fromfunction

	def __str__(self):
		if (self.fromfunction is None):
			return "%s %s\n" % (self.time,self.msg)
		else:
			return "%s %s (function: %s)\n" % (self.time,self.msg,self.fromfunction)

class Log(object):
	def __init__(self,printit=True):
		self.info = []
		self._printit = printit

	def __call__(self,msg=None,func=None,level=None):
		if (msg is None):
			print self.tostr()
			return
		if (level is None or level>=g_loglevel):
			msg = Message(msg,func)
			if (self._printit): print str(msg).strip()
			self.info.append(msg)

	def clean(self):
		del self.info
		self.info = []

	def disablePrintScreen(self):
		self._printit = False

	def enablePrintScreen(self):
		self._printit = True

	def tofile(self,fname):
		f=open(fname,"a")
		f.write("## Writing new log\n")
		f.write("%s" % self.tostr())
		f.close()

	def tostr(self):
		s = ""
		for i in self.info:
			if (i.fromfunction is not None):
				s += "%s %s (function: %s)\n" % (i.time,i.msg,i.fromfunction)
			else:
				s += "%s %s\n" % (i.time,i.msg)
		return s
		

logbook=Log()
logbook("Starting new log")
