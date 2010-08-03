from distutils.core import setup, Extension, Command
from distutils.command.install import install,INSTALL_SCHEMES
from distutils.command.build import build
from distutils.util import get_platform
from distutils import log as dlog
from subprocess import Popen
from gimmemotifs.utils import which
from gimmemotifs.tools import *
from glob import glob
import os
import sys
import shutil
from stat import ST_MODE
import inspect

CONFIG_NAME = "gimmemotifs.cfg" 
DESCRIPTION  = """GimmeMotifs is a motif prediction pipeline. 
"""

DEFAULT_PARAMS = {
	"analysis": "medium",
	"fraction": 0.2,
	"abs_max": 1000,
	"width": 200,
	"lwidth": 500,
	"pvalue": 0.001,
	"enrichment": 1.5,
	"background": "genomic_matched,random",
	"genome": "hg18",
	"tools": "MDmodule,Weeder,MotifSampler",
	"available_tools": "Weeder,MDmodule,MotifSampler,gadem,meme,trawler,WannaMotif,Improbizer,MoAn,BioProspector",
	"cluster_threshold": "0.95",
	"use_strand": False
}

MOTIF_CLASSES = ["MDmodule", "Meme", "Weeder", "Gadem", "MotifSampler", "Trawler", "Improbizer", "MoAn", "BioProspector"]

# Included binaries after compile
MOTIF_BINS = {
	"Meme": "src/meme_4.4.0/src/meme.bin",
	"MDmodule": "src/MDmodule/MDmodule",
	"BioProspector": "src/BioProspector/BioProspector",
	"MoAn": "src/MoAn/moan",
	"Gadem": "src/GADEM_v1.3/src/gadem"
}

data_files=[
	('gimmemotifs/templates', ['templates/cluster_template.kid', 'templates/report_template.kid', 'templates/report_template_v2.kid', 'templates/cluster_template_v2.kid']),
	('gimmemotifs/score_dists', ['score_dists/total_wic_mean_score_dist.txt']),
	('gimmemotifs/genes', ['genes/hg18.bed', 'genes/mm9.bed', 'genes/xenTro2.bed']),
	('gimmemotifs/bg', ['bg/hg18.MotifSampler.bg', 'bg/mm9.MotifSampler.bg', 'bg/xenTro2.MotifSampler.bg']),
	('gimmemotifs/motif_databases', ['motif_databases/jaspar.pfm']),
	('gimmemotifs/genome_index', ['genome_index/README.txt'])
]


# Fix for install_data, add share to prefix (borrowed from Dan Christiansen) 


for platform, scheme in INSTALL_SCHEMES.iteritems():
	if platform.startswith('unix_'):
		if scheme['data'][0] == '$' and '/' not in scheme['data']:
			scheme['data'] = os.path.join(scheme['data'], 'share')

dlog.info("checking dependencies")
try:
	import pp
	import matplotlib
	import kid
	import scipy
	import numpy
except ImportError, inst:
	print "Error: required dependency not found!"
	print inst
	sys.exit()	

class build_tools(Command):
	description = "compile all included motif prediction tools"

	def initialize_options(self):
		self.build_base = None
		self.plat_name = None
		self.build_tools_dir = None

	def finalize_options(self):	
		if self.plat_name is None:
			self.plat_name = get_platform()
		self.set_undefined_options('build',('build_base', 'build_base'))
		plat_specifier = ".%s-%s" % (self.plat_name, sys.version[0:3])
		self.build_tools_dir = os.path.join(self.build_base, 'tools' + plat_specifier)
		self.set_undefined_options('build',('custom_build', 'build_tools_dir'))
	
	def run(self):
		from compile_externals import compile_all
		if not os.path.exists(self.build_tools_dir):
			os.mkdir(self.build_tools_dir)

		# Try to compile everything
		compile_all()

		# Copy everything that has been compiled
		for bin in MOTIF_BINS.values():
			if os.path.exists(bin):
				shutil.copy(bin, self.build_tools_dir)

		# Copy seqlogo
		if os.path.exists("src/weblogo"):
			dlog.info("building seqlogo")
			patterns = ["src/weblogo/logo.*", "src/weblogo/template.*", "src/weblogo/seqlogo"]
			for p in patterns:
				for file in glob(p):
					shutil.copy(file, self.build_tools_dir)

		# Copy trawler
		if os.path.exists("src/trawler_standalone-1.2"):
			dlog.info("building trawler")
			if os.path.exists(os.path.join(self.build_tools_dir, "trawler")):
				shutil.rmtree(os.path.join(self.build_tools_dir, "trawler"))
			shutil.copytree("src/trawler_standalone-1.2", os.path.join(self.build_tools_dir, "trawler"))

class build_config(Command):
	description = "create a rudimentary config file"
	
	def initialize_options(self):
		self.build_cfg = None
		self.build_base = None
		self.build_tools_dir = None
	
	def finalize_options(self):	
		self.set_undefined_options('build', ('build_base', 'build_base'))
		self.set_undefined_options('build_tools', ('build_tools_dir', 'build_tools_dir'))
		#self.set_undefined_options('install', ('install_data', 'install_dir'))
		self.build_cfg = os.path.join(self.build_base, "cfg")

	def run(self):
		if not os.path.exists(self.build_cfg):
			os.mkdir(self.build_cfg)

		from gimmemotifs.config import MotifConfig
		cfg = MotifConfig(use_config="cfg/gimmemotifs.cfg.example")
		
		dlog.info("locating motif programs")
		available = []
		for program in MOTIF_CLASSES:
			# Get class
			m = eval(program)()
			cmd = m.cmd
			
			### ugly, fixme :)
			if cmd == "trawler.pl":
				cmd = "trawler/bin/trawler.pl"

			bin = ""
			if os.path.exists(os.path.join(self.build_tools_dir, cmd)):
				bin = os.path.join(self.build_tools_dir, cmd)
				dlog.info("using included version of %s: %s" % (program, bin))
			else:
				### ugly, fixme :)
				if cmd == "trawler/bin/trawler.pl":
					cmd = "trawler.pl"
				if program in MOTIF_BINS.keys():
					dlog.info("could not find compiled version of %s" % program)
				bin = which(cmd)
				if bin:
					dlog.info("using installed version of %s: %s" % (program, bin))
				else:
					dlog.info("not found: %s" % program)
			
			### Some more ugly stuff
			if bin:
				dir = bin.replace(m.cmd,"")
				if program == "Weeder":
					dir = bin.replace("weederTFBS.out","")
				elif program == "Meme":
					dir = bin.replace("bin/meme.bin", "").replace("meme.bin", "")
				elif program == "Trawler":
					dir = bin.replace("bin/trawler.pl", "")

				available.append(m.name)
				cfg.set_program(m.name, {"bin":bin, "dir":dir})

		# Weblogo
		bin = ""
		seq_included = os.path.join(self.build_tools_dir, "seqlogo")
		if os.path.exists(seq_included):
			bin = seq_included
			dlog.info("using included version of weblogo: %s" % seq_included)
		else:
			bin = which("seqlogo")
			dlog.info("using installed version of seqlogo: %s" % (bin))
		if bin:
			cfg.set_seqlogo(bin)
		else:
			dlog.info("couldn't find seqlogo")
		
		# Set the available tools in the config file
		DEFAULT_PARAMS["available_tools"] = ",".join(available)
		DEFAULT_PARAMS["tools"] = ",".join(available)
		cfg.set_default_params(DEFAULT_PARAMS)

		# Write (temporary) config file
		config_file = os.path.join(self.build_cfg, "%s" % CONFIG_NAME)
		dlog.info("writing (temporary) configuration file: %s" % config_file)
		f = open(config_file, "wb")
		cfg.write(f)
		f.close()

	def get_outputs(self):
		return self.outfiles or []

class install_tools(Command):
	description = "install (compiled) motif prediction tools"
	
	def initialize_options(self):
		self.tools_dir = None
		self.install_dir = None
		self.install_tools_dir = None
	
	def finalize_options(self):	
		self.set_undefined_options('build_tools', ('build_tools_dir', 'tools_dir'))
		self.set_undefined_options('install', ('install_data', 'install_dir'))
		self.install_tools_dir = os.path.join(self.install_dir, "gimmemotifs/tools")

	def run(self):
		dir = "src/Algorithm-Cluster-1.49/"
		if os.path.exists(os.path.join(dir, "Makefile")):
			Popen(["make","install"], cwd=dir, stdout=PIPE).communicate()

		dst = os.path.join(self.install_dir, "gimmemotifs/tools")
		self.outfiles = self.copy_tree(self.tools_dir, self.install_tools_dir)
	def get_outputs(self):
		return self.outfiles or []

class install_config(Command):
	description = "create and install a customized configuration file"

	def remove_nonsense(self, dir):
		if dir.find("BUILDROOT") != -1:
			components = os.path.normpath(os.path.abspath(dir)).split(os.sep)
			for i in range(len(components)):
				if components[i] == "BUILDROOT":
					return os.path.sep.join([""] + components[i + 2:])
		elif dir.find("debian") != -1:
			components = os.path.normpath(os.path.abspath(dir)).split(os.sep)
			for i in range(len(components)):
				if components[i] == "debian":
					return self.remove_nonsense(os.path.sep.join([""] + components[i + 2:]))
			
		return dir


	def initialize_options(self):
		self.build_base = None
		self.install_dir = None
		self.build_cfg = None
		self.build_tools_dir = None
		self.install_tools_dir = None

	def finalize_options(self):
		self.set_undefined_options('build', ('build_base', 'build_base'))
		self.set_undefined_options('install', ('install_data', 'install_dir'))
		self.set_undefined_options('build_config', ('build_cfg', 'install_dir'))
		self.set_undefined_options('build_tools', ('build_tools_dir', 'build_tools_dir'))
		self.set_undefined_options('install_tools', ('install_tools_dir', 'install_tools_dir'))
	
	def run(self):
		from gimmemotifs.config import MotifConfig
		
		cfg = MotifConfig(use_config=self.build_cfg)

		data_dir = self.remove_nonsense(os.path.abspath(self.install_dir))
		
		cfg.set_template_dir(os.path.join(data_dir, 'gimmemotifs/templates'))
		cfg.set_gene_dir(os.path.join(data_dir, 'gimmemotifs/genes'))
		cfg.set_score_dir(os.path.join(data_dir, 'gimmemotifs/score_dists'))
		cfg.set_index_dir(os.path.join(data_dir, 'gimmemotifs/genome_index'))
		cfg.set_motif_dir(os.path.join(data_dir, 'gimmemotifs/motif_databases'))
		cfg.set_bg_dir(os.path.join(data_dir, 'gimmemotifs/bg'))
		cfg.set_tools_dir(os.path.join(data_dir, 'gimmemotifs/tools'))
		
		final_tools_dir = self.remove_nonsense(self.install_tools_dir)
		for program in MOTIF_CLASSES:
			m = eval(program)()
			if cfg.is_configured(m.name):
				bin = cfg.bin(m.name).replace(self.build_tools_dir, final_tools_dir) 
				dir = cfg.dir(m.name)
				if dir:
					dir = dir.replace(self.build_tools_dir, final_tools_dir)
				cfg.set_program(m.name, {"bin":bin, "dir":dir})
			
		dir = cfg.get_seqlogo()
		dir = dir.replace(self.build_tools_dir, final_tools_dir)
		cfg.set_seqlogo(dir)

		# Use a user-specific configfile if any other installation scheme is used
#		if os.path.abspath(self.install_dir) == "/usr/share":
		config_file = os.path.join(self.install_dir, "gimmemotifs/%s" % CONFIG_NAME)
		self.outfiles = [config_file] 
#		else:
#			config_file = os.path.expanduser("~/.%s" % CONFIG_NAME)
		
		if os.path.exists(config_file):
			new_config = config_file + ".tmp"
			dlog.info("INFO: Configfile %s already exists!\n      Will create %s, which contains the new config.\n      If you want to use the newly generated config you can move %s to %s, otherwise you can delete %s.\n" % (config_file, new_config, new_config, config_file, new_config))

			f =  open(new_config, "wb")
			cfg.write(f)
		else: 
			dlog.info("writing configuration file %s" % config_file)
			f =  open(config_file, "wb")
			cfg.write(f)
		
		dlog.info("edit %s to further configure GimmeMotifs." % config_file)
	
	def get_outputs(self):
		return self.outfiles or []

class custom_build(build):
	def run(self):
		build.run(self)
		self.run_command('build_tools')
		self.run_command('build_config')

class custom_install(install):
	sub_commands = install.sub_commands + [
			('install_tools', lambda self: True),
			('install_config', lambda self: True)
			]

	# Make sure we install in the correct locations on Ubuntu
	def finalize_options(self):
		install.finalize_options(self)
		if self.install_data == "/usr":
			self.install_data = "/usr/share"
		if self.install_data.endswith("/usr"):
			parts = self.install_data.split(os.sep)
			if parts[-3] == "debian":
				self.install_data = os.path.join(self.install_data, "share")

	
	def run(self):
		install.run(self)
	
module1 = Extension('gimmemotifs.c_metrics', sources = ['gimmemotifs/c_metrics.c'], libraries = ['gsl', 'gslcblas'])

setup (name = 'gimmemotifs',
		cmdclass={"build":custom_build, 
							"build_tools":build_tools,
							"build_config":build_config,
							"install":custom_install, 
							"install_tools":install_tools,
							"install_config":install_config,
							},
		version = '0.60',
		description = DESCRIPTION,
		author='Simon van Heeringen',
		author_email='s.vanheeringen@ncmls.ru.nl',
		url='http://www.ncmls.eu/bioinfo/gimmemotifs/',
		license='MIT',
		packages=['gimmemotifs'],
		ext_modules = [module1],
		scripts=[
			'scripts/add_organism.py',
			'scripts/generate_sequences.py',
			'scripts/closest_motif_match.py',
			'scripts/motif_cluster.py',
			'scripts/create_genome_index.py',
			'scripts/gimme_motifs.py',
			'scripts/motif_roc.py',
			'scripts/motif_roc_metrics.py',
			'scripts/motif_localization_plots.py',
			'scripts/pwm2logo.py',
			'scripts/track2fasta.py',
			'scripts/pwmscan.py',
			],
		data_files=data_files
)
