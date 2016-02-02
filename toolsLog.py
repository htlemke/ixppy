from __future__ import print_function
import utilities
import sys

g_loglevel=1

class Message(object):
  def __init__(self,msg,fromfunction=None,showTimestamp=True):
    self.msg=msg
    self.time = utilities.now()
    self.fromfunction=fromfunction
    self.showTimestamp=showTimestamp

  def __str__(self):
    if self.showTimestamp:
      msg = "%s " % self.time
    else:
      msg = ""
    if (self.fromfunction is None):
      msg +="%s" % self.msg
    else:
      msg += "%s (function: %s)" % (self.msg,self.fromfunction)
    return msg

class Log(object):
  def __init__(self,printit=True):
    self.info = []
    self._printit = printit

#  def __call__(self,*args,func=None,level=None,end="\n",printit=None):
  def __call__(self,*messages,**kw):
    func  = kw.get("func",None)
    level = kw.get("level",None)
    end   = kw.get("end","\n")
    showTimestamp = kw.get("time",True)
    printit = kw.get("printit",None)
    if printit is None: printit = self._printit
    if (len(messages)==0):
      print(self.tostr()+end)
      return
    else:
      if (level is None or level>=g_loglevel):
        msg = ["%s" % m for m in messages]
        msg = " ".join(msg)
        msg = Message(msg,func,showTimestamp=showTimestamp)
        self.info.append(msg)
        if (self._printit):
          print(msg,end=end)
          sys.stdout.flush()

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
    


def msg(s,newline=True):
  sys.stdout.write(s)
  if (newline):
    sys.stdout.write("\n")
  sys.stdout.flush()

class codeBlock(object):
  def __init__(self,what,level=0):
    self.mem = mem()
    self.t0=time.time()
    self.space = "  "*level
    s=self.space+what+" (%.1f Gb)..." % (self.mfree())
    msg(s,newline=False)
  def done(self):
    t=time.time()-self.t0
    msg("...  done (%.1f Gb, %.1f sec)" % (self.mfree(),t))
  def mfree(self):
    
    return self.mem.free/1024./1024./1024.

logbook=Log()
logbook("Starting new log")
