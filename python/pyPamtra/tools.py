# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import traceback
import warnings
import sys

try: 
  import paramiko 
except: 
  warnings.warn("paramiko not available", Warning)


class sftp2Cluster(object):

  def __init__(self,machinename, username):
    self.ssh = paramiko.SSHClient()
    self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    self.ssh.connect(machinename, username=username)
    self.sftp = self.ssh.open_sftp()
    
  def __del__(self):
    self.ssh.close()
  
  def put(self, filename, data):
    try:
        self.sftp.mkdir(os.path.dirname(filename))
    except IOError:
        pass
    f = self.sftp.open(filename, 'w')
    f.write(data)
    f.close()

  def rm(self, filename,silent=True):
    try: 
      self.sftp.remove(filename)
      print "removed %s"%(filename)
    except IOError:
      if silent:
        pass
      else:
        raise IOError

  def mv(self, oldFile, newFile):
    self.sftp.rename(oldFile, newFile)

def formatExceptionInfo(maxTBlevel=5):
  cla, exc, trbk = sys.exc_info()
  excName = cla.__name__
  try:
    excArgs = exc.__dict__["args"]
  except KeyError:
    excArgs = "<no args>"
  excTb = traceback.format_tb(trbk, maxTBlevel)
  return (excName, excArgs, excTb)

  