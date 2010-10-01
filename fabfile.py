# Borrows heavily from code by whiteinge http://gist.github.com/234001
from __future__ import with_statement
from fabric.api import *
from fabric.contrib.console import confirm
import subprocess
import time
import socket

env.hosts = ['localhost:2223']

def _local(*args, **kwargs):
	"""A wrapper around fabric.local() that takes a list. It will take care of
	the shell quoting for you. This makes interspersing variables a bit easier.
	
	E.g.: instead of::

		local('echo "Hello, this has spaces in it.")

	Use::
											        
		_local(['echo', 'This has spaces in it.'])

	"""
	return local(subprocess.list2cmdline(*args), **kwargs)


def _vbox_state(machine):
	"""Get the state of the specified machine."""
	vminfo = subprocess.Popen(['VBoxManage', '-q', 'showvminfo', machine, 
		'--machinereadable'], stdout=subprocess.PIPE)
	result = filter(lambda x: 'VMState="' in x, vminfo.stdout.readlines()).pop()
	return result.split('"')[1]

def _vbox_start(machine, port=2222):
	""" Start the virtual machine and wait until it's up and running."""
	
	# Check if the machine is already running
	if _vbox_state(machine) == "running":
		return 

	# Start the virtual machine
	_local(['VBoxManage','-q','startvm',machine,'--type','headless'])

	while True:
		s = None
		try:
			s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
			s.settimeout(3.0)
			s.connect(('localhost', port))
			response = s.recv(1024)
		except (socket.error, socket.timeout), e:
			pass
		else:
			if 'SSH' in response:
				break

			time.sleep(3)
		finally:
			if s:
				s.close()

	
def _vbox_stop(machine):
	_local(['VBoxManage','-q','controlvm',machine,'poweroff'])

	while True:
		state = _vbox_state(machine)

		if state in ['stopped', 'poweroff', 'saved']:
			break
		
		print state
		time.sleep(3)

def _vbox_load_snapshot(machine, name, force=True):
	""" Load a specific snapshot for a machine"""
	
	if _vbox_state(machine) == "running":
		if force:
			_vbox_stop(machine)
		else:
			raise "Can't restore state of running machine!"

	_local(['VBoxManage','-q','snapshot',machine,'restore', name])

def deploy():
	#put('/tmp/my_project.tgz', '/tmp/')
	with cd('git/gimmemotifs'):
		run('git pull origin master')
		#run('python make_deb.py')
	get('git/gimmemotifs/build/debian/gimmemotifs_0.61-1_i386.deb', '/tmp/')

def build():
	# Load snapshot
	# Start machine
	snapshot = "autobuild2"
	machine = "Ubuntu-9.10-x86_64"
	
	_vbox_load_snapshot(machine, snapshot)
	_vbox_start(machine, port=2223)
	
	with cd('git/gimmemotifs'):
		run('git pull origin master')
		run('python make_deb.py')
	get('git/gimmemotifs/build/debian/gimmemotifs/gimmemotifs_0.61-1_x86_64.deb', '/tmp/')
	
	_vbox_stop(machine)

if __name__ == "__main__":
	#snapshot = "c0bf134f-8426-4e2b-bd2e-f6857c302e52"
	#machine = "Ubuntu-9.10-i386"
	#snapshot = "db8acfda-1bd5-48e7-a2b6-ce73059d6e16"
	snapshot = "autobuild"
	machine = "Ubuntu-9.10-x86_64"
	
	#_vbox_stop(machine)

