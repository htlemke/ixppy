import sys
import datetime
import subprocess
import os
from toolsVarious import iterfy

preamble  = 'import sys\nimport ixppy\nfrom ixppy.lsfHelper import import_file\n'
queue = 'psnehq'

def import_file(full_path_to_module,glb):
  #try:
        import os
        module_dir, module_file = os.path.split(full_path_to_module)
        module_name, module_ext = os.path.splitext(module_file)
        save_cwd = os.getcwd()
        os.chdir(module_dir)
	os.system('pwd')
	print module_name
	module_obj = __import__(module_name)
	#exec('import '+module_name)
	module_obj.__file__ = full_path_to_module
        glb[module_name] = module_obj
	#raise NotImplementedError('Use the source, luke!')
        os.chdir(save_cwd)
	return module_name
    #except Exception,e:
	#print e
	#raise ImportError

def get_module_name(full_path_to_module):
    try:
        import os
        module_dir, module_file = os.path.split(full_path_to_module)
        module_name, module_ext = os.path.splitext(module_file)
	return module_name
    except Exception,e:
        print e
        raise ImportError

def writePyFile(funcs,experiment=None,save=True,identifier=''):
  timestr = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
  filename = 'lsf_autocall_'+timestr+'.py'
  filestr = preamble
  if experiment==None:
    filestr += 'experiment = sys.argv[1]\nprint "got experiment %s"%experiment\n'
    filestr += 'runnos = [int(tr) for tr in sys.argv[2:]]\nprint "got run %s"%runnos\n'
  else:
    filestr += 'runnos = [int(tr) for tr in sys.argv[1:]]\nprint "got run %s"%runnos\n'
    #filestr += 'runno = int(sys.argv[1])\nprint "got run %d"%runno\n'
  modnames = []
  for func in funcs:
    modFina = sys.modules[func.__module__].__file__
    modFina = os.path.abspath(modFina)
    modnames.append(get_module_name(modFina))
    filestr += 'import_file(\'%s\',globals())\n'%modFina
  if experiment==None:
    filestr += 'cachefina = experiment+\'run\'+\'-\'.join(str(run) for run in runnos)+\'_%s_lsfHelper.ixp.h5\'\n' %identifier
    filestr += 'd = eval(\'ixppy.dataset((\\\'%s\\\',%d),ixpFile=%s)\'%(experiment,runno,cachefina))\n'
  else:
    filestr += 'cachefina = \'%s\'+\'run\'+\'-\'.join(str(run) for run in runnos)+\'_%s_lsfHelper.ixp.h5\'\n'%(experiment,identifier)
    filestr += 'd = ixppy.dataset((\'%s\',runnos),ixpFile=cachefina)\n'%experiment
  for func,modname in zip(funcs,modnames):
    funcname = func.__name__
    filestr += modname+'.'+funcname+'(d)\n'
    if save:
      filestr += 'd.save(force=True)\n'
    filestr += 'print \'Batch job completed successfully.\'\n'
  h = file(filename,'w')
  h.write(filestr)
  h.close()
  return filename

def applyFileOnRun(fina,runno,experiment=None):
  execstr = 'bsub -q ' + queue + ' '
  execstr+= '-o ' + os.path.splitext(fina)[0]
  if not experiment==None:
    execstr+= '_'+experiment
  for tr in runno:
    execstr+= '_run%04d'%tr
  execstr+= '.out '
  execstr+=' python '
  execstr+=fina+' '
  if not experiment==None:
    execstr+= experiment+' '
  for run in runno:
    execstr+= str(run) + ' '
  result = subprocess.check_output(execstr, shell=True)
  #os.system(execstr)
  return result

def applyOnRunlist(input,runnos,experiment,save=True,identifier=''):
  if type(input) is str:
    fina = input
  elif type(input) is list:
    funcs = input
    fina = writePyFile(funcs,experiment=experiment,save=save,identifier=identifier)
  returns = []
  for runno in runnos:
    if len(iterfy(runno))>1:
      returns.append(applyFileOnRun(fina,runno))
    else:
      returns.append(applyFileOnRun(fina,runno))
  return returns





  
  
