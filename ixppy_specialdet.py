import numpy
import h5py
import os
from socket import gethostname
import dateutil
import sys
import numpy as np
import ixppy

class timeTool(object):
  def __init__(self,config):
    self.config = config
    self._dataset_data = [['Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:AMPL/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:AMPLNXT/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:FLTPOS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:FLTPOSFWHM/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:FLTPOS_PS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.2/TTSPEC:REFAMPL/data',],
                          ['Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:AMPL/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:AMPLNXT/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:FLTPOS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:FLTPOSFWHM/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:FLTPOS_PS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.1/TTSPEC:REFAMPL/data',],
                          ['Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:AMPL/data',
                          'Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:AMPLNXT/data',
                          'Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:FLTPOS/data',
                          'Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:FLTPOSFWHM/data',
                          'Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:FLTPOS_PS/data',
                          'Epics::EpicsPv/EpicsArch.0:NoDevice.0/TTSPEC:REFAMPL/data',],
                          
                          ['Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:AMPL/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:AMPLNXT/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:FLTPOS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:FLTPOSFWHM/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:FLTPOS_PS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TTSPEC:REFAMPL/data',]
                          ,['Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:AMPL/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:AMPLNXT/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:FLTPOS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:FLTPOSFWHM/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:FLTPOS_PS/data',
                          'Epics::EpicsPv/XppEndstation.0:Opal1000.0/TIMETOOL:RAWPOS/data',]]
    self._dataset_fieldnames = ['ampl','amplnxt','fltpos','fltposfwhm'
                                ,'fltpos_ps','rawpos_or_refampl']
    
  def rdAllData(self):
    self.config.base._rdScanPar()
    data = []

    # find data format
    for datasets in self._dataset_data:
      try:
        dsetstring = self.config.base._mkScanStepDataSetString(datasets[0],0)
        dset = self.config.fileHandle[dsetstring]
        tdatasets = datasets
        print "found good timetool dataset %s" %datasets[0]
        break
      except:
        continue

    for stepNo in range(len(self.config.base._controlPv)):
      dsetNO = 0
      fieldnames = ''
      tStepData = []
      #for datasets in self._dataset_data:
        #try:
      for tdetdset in tdatasets:
        try:
          dsetstring = self.config.base._mkScanStepDataSetString(tdetdset,stepNo)
          dset = self.config.fileHandle[dsetstring]
          tempdat = ixppy.rdHdf5dataFromDataSet(dset)
          tStepData.append(tempdat['value'])
          fieldnames+= self._dataset_fieldnames[dsetNO] + ', '
          dsetNO+=1
        except:
          dsetNO+=1
          print "sht"
          continue
        #except:
          #print "NB: older version of Timing too data format!"
          #continue
      fieldnames = fieldnames[:-2]
      #dtypes = [np.dtype(xx[0]) for xx in tStepData]
      data.append(np.core.records.fromarrays(np.array(tStepData),names=fieldnames))
      #data.append(np.array(tStepData),names=fieldnames))

    return data



  def rdStepData(self):
    print 'Not implemented!!'


class example(object):
  def __init__(self,config):
    self.config = config
  def rdDetAllData(self):
    return data

  def rdStepData(self):
    return data


class eventCode(object):
  def __init__(self,config):
    self.config = config
    self._dataset_data =  []
    self._dataset_data.append('EvrData::DataV3/NoDetector.0:Evr.0/evrData')
    self._dataset_data.append('EvrData::DataV3/NoDetector.0:Evr.0/data')
    self._dataset_config = []
    self._dataset_config.append('/Configure:0000/EvrData::ConfigV7/NoDetector.0:Evr.0/eventcodes')
    self._dataset_config.append('/Configure:0000/EvrData::ConfigV5/NoDetector.0:Evr.0/eventcodes')

  def rdAllData(self):
    #try:
      self.config.base._rdScanPar()
      data = []
      for conf_dsetname in self._dataset_config:
        try:
          config_dset = self.config.fileHandle[conf_dsetname]
          break
        except:
          pass

      codes = config_dset['code']
      names = ['code%d'%code for code in codes]
      for stepNo in range(len(self.config.base._controlPv)):
        dsetNO = 0
        fieldnames = ''
        tStepData = []
        for dataset_string in self._dataset_data:
          try:
            dsetstring = self.config.base._mkScanStepDataSetString(dataset_string,stepNo)
            dset = self.config.fileHandle[dsetstring]
            break
          except:
            pass

        tdat = dset.value
        boolmat = np.zeros([len(codes),len(tdat)],dtype=bool)
        for shot,shotNO in zip(tdat,range(len(tdat))):
          for code,codeNO in zip(codes,range(len(codes))):
            if code in shot[0][0]:
              boolmat[codeNO,shotNO] = True

        data.append(np.core.records.fromarrays(boolmat,names=names))

      return data
    #except:
      #print "Could not eread event code data, could be relatet to the datqa format, evr data can presently only be read by a patched version of h5py on lcls servers."


  
  #def rdDetAllData(self):
    #return data

  def rdStepData(self):
    return data
