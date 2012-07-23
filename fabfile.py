# Script to compile GimmeMotifs for packaging for several different distributions
# A VirtualBox image needs to be present for every "host"
# Borrows some code from whiteinge http://gist.github.com/234001
#
from __future__ import with_statement
from fabric.api import *
from fabric.contrib.console import confirm
from fabric import state
import subprocess
import time
import socket
import os


DEFAULT_PORT = 22

HOST_DATA = {
	'192.168.56.10': ["Fedora 16 32bit", 
										"autobuild",
										"python make_rpm.py",
										"rpmbuild/RPMS/i686/gimmemotifs-%s-1.i686.rpm",
										"gimmemotifs-%s-1.i686.rpm"],
	'192.168.56.11': ["Fedora 16 64bit", 
										"autobuild",
										"python make_rpm.py 2>&1",
										"rpmbuild/RPMS/x86_64/gimmemotifs-%s-1.x86_64.rpm",
										"gimmemotifs-%s-1.x86_64.rpm"],
	'192.168.56.12': ["Debian 6.0 32bit", 
										"autobuild",
										"python make_deb.py 2>&1",
										"git/gimmemotifs/build/debian/gimmemotifs_%s_i386.deb",
										"gimmemotifs_%s_i386.debian.deb"],
	'192.168.56.13': ["Debian 6.0 64bit", 
										"autobuild",
										"python make_deb.py 2>&1",
										"git/gimmemotifs/build/debian/gimmemotifs_%s_amd64.deb",
										"gimmemotifs_%s_amd64.debian.deb"],
	'192.168.56.14': ["Ubuntu 11.10 32-bit", 
										"autobuild",
										"python make_deb.py 2>&1",
										"git/gimmemotifs/build/debian/gimmemotifs_%s_i386.deb",
										"gimmemotifs_%s_i386.ubuntu.deb"],
	'192.168.56.15': ["Ubuntu 11.10 64bit", 
										"autobuild",
										"python make_deb.py 2>&1",
										"git/gimmemotifs/build/debian/gimmemotifs_%s_amd64.deb",
										"gimmemotifs_%s_amd64.ubuntu.deb"]
}

env.hosts = HOST_DATA.keys()

def _local(*args, **kwargs):
	"""A wrapper around fabric.local() that takes a list. It will take care of
	the shell quoting for you. This makes interspersing variables a bit easier.
	"""
	return local(subprocess.list2cmdline(*args), **kwargs)

def _vbox_state(machine):
	"""Get the state of the specified machine."""
	vminfo = subprocess.Popen(['VBoxManage', '-q', 'showvminfo', machine, 
		'--machinereadable'], stdout=subprocess.PIPE)
	result = filter(lambda x: 'VMState="' in x, vminfo.stdout.readlines()).pop()
	return result.split('"')[1]

def _vbox_start(machine, server="localhost", port=2222):
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
			s.connect((server, port))
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
	""" Stop a virtual macbhine"""

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

def _get_version():
	""" Get package version, based on local repository """
	version  = subprocess.Popen(["python","setup.py","--version"], stdout=subprocess.PIPE).communicate()[0].strip()
	return version
		
def _build_general(machine, snapshot, server, port, build_cmd, build_file, local_file):
	""" The scaffold of the build sequence. """
	# Load snapshot
	_vbox_load_snapshot(machine, snapshot)
	# Start machine
	_vbox_start(machine, server=server, port=port)
	# Wait to make sure networking is up and running
	time.sleep(5)
	
	# Checkout latest version and compile
	with cd('git/'):
		run('git clone git://github.com/simonvh/gimmemotifs.git')
	with cd('git/gimmemotifs'):
		run(build_cmd, pty=False)
	
	# Download the package
	file = os.path.join("dist", local_file)
	get(build_file, file)

	# Rename if fabric insists on including the hostname in the file
	if os.path.exists("%s.%s" % (file, state.env.host)):
		os.rename("%s.%s" % (file, state.env.host), file)
	
	# And we're done with this virtual machine
	_vbox_stop(machine)

def _test_ubuntu(machine, snapshot, server, port):
	""" Test GimmeMotifs installation on Ubuntu """
	
	# Load snapshot
	_vbox_load_snapshot(machine, snapshot)
	# Start machine
	_vbox_start(machine, server=server, port=port)
	# Wait to make sure networking is up and running
	time.sleep(5)
	
	# Checkout latest version and compile
	#with cd('Downloads/'):
	#	run('git clone git://github.com/simonvh/gimmemotifs.git')
	#with cd('git/gimmemotifs'):
	#	run(build_cmd)
	
	# And we're done with this virtual machine
	_vbox_stop(machine)



def build_packages():
	""" Build packages for all specified distributions"""

	# Get current host
	hostname = state.env.host
	
	# Get parameters for this machine / host
	machine, snapshot, build_cmd, build_file, local_file = HOST_DATA[hostname]
	
	# Get the current version
	version = _get_version()
	# Replace %s with current version
	build_file = build_file % version
	local_file = local_file % version
	
	# Build the package on each host
	_build_general(machine, snapshot, hostname, DEFAULT_PORT, build_cmd, build_file, local_file)
