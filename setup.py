from distutils.core import setup, Extension
from distutils.command.install import INSTALL_SCHEMES
from distutils.command.build import build
from  distutils.command.install_data import install_data
#from distutils import log
from gimmemotifs.utils import which
from gimmemotifs.tools import *
from glob import glob
import os
import sys
from shutil import copyfile
from stat import ST_MODE

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


print "Checking dependencies"
try:
	import pp
	import matplotlib
	import kid
	import scipy
	import numpy
except ImportError as inst:
	print "Error: required dependency not found!"
	print inst
	sys.exit()	

class custom_build(build):
	def run(self):
		build.run(self)
		self.compile_external()
	
	def compile_external(self):
		from compile_externals import compile_all
		compile_all()

class install_plus_config(install_data):
	def run(self):
		#print help(self)
		#print self.distribution.data_files
		#print self.install_data
		#sys.exit()
		self.write_config()
		install_data.run(self)
		self.chmod_tools()
	
	def chmod_tools(self):
		data_dir = os.path.abspath(self.install_dir)
		dir = os.path.join(data_dir, "gimmemotifs/tools")
		for file in os.listdir(dir):
			mode = ((os.stat(os.path.join(dir,file))[ST_MODE]) | 0555) & 07777		
			#log.info("changing mode of %s to %o", file, mode)
			print os.path.join(dir,file)

			os.chmod(os.path.join(dir, file), mode)


	def write_config(self):
		from gimmemotifs.config import MotifConfig
		cfg = MotifConfig(use_config="cfg/gimmemotifs.cfg.example")

		data_dir = os.path.abspath(self.install_dir)
		cfg.set_template_dir(os.path.join(data_dir, 'gimmemotifs/templates'))
		cfg.set_gene_dir(os.path.join(data_dir, 'gimmemotifs/genes'))
		cfg.set_score_dir(os.path.join(data_dir, 'gimmemotifs/score_dists'))
		cfg.set_index_dir(os.path.join(data_dir, 'gimmemotifs/genome_index'))
		cfg.set_motif_dir(os.path.join(data_dir, 'gimmemotifs/motif_databases'))
		cfg.set_bg_dir(os.path.join(data_dir, 'gimmemotifs/bg'))
		cfg.set_tools_dir(os.path.join(data_dir, 'gimmemotifs/tools'))
		
		print
		print "Locating motif programs"
		included_bins = []
		available = []
		for program in MOTIF_CLASSES:
			m = eval(program)()
			cmd = m.cmd
			bin = which(cmd)
			
			if bin:
				print "Found %s in %s" % (m.name, bin)
			else:
				if program in MOTIF_BINS.keys() and os.path.exists(MOTIF_BINS[program]):
					print "Using included compiled version of %s" % program
					bin = MOTIF_BINS[program]
					copyfile(bin, os.path.join("tools", cmd))
					included_bins.append(os.path.join("tools", cmd))
					bin = os.path.join(data_dir, "gimmemotifs/tools", cmd)
				else:
					print "%s not installed!" % program
			
			if bin:
				available.append(m.name)
				dir = None
				if program == "Weeder":
					dir = bin.replace("weederTFBS.out","")
				elif program == "Meme":
					dir = bin.replace("bin/meme.bin", "")
				elif program == "Trawler":
					dir = bin.replace("bin/trawler.pl", "")
		
				cfg.set_program(m.name, {"bin":bin, "dir":dir})
			else:
				print "Couldn't find %s" % m.name
	
		if included_bins:
			self.distribution.data_files.append	(("gimmemotifs/tools", included_bins))
		
		print
		print "Trying to locate seqlogo"
		bin = which("seqlogo")
		if bin:
			print "Found seqlogo in %s" % (bin)
			cfg.set_seqlogo(bin)
		else:
			print "Couldn't find seqlogo"
		print
		
		DEFAULT_PARAMS["available_tools"] = ",".join(available)
		DEFAULT_PARAMS["tools"] = ",".join(available)
		cfg.set_default_params(DEFAULT_PARAMS)
		
		# Use a user-specific configfile if any other installation scheme is used
		if os.path.abspath(self.install_dir) == "/usr/share":
			config_file = "/usr/share/gimmemotifs/%s" % CONFIG_NAME
		else:
			config_file = os.path.expanduser("~/.%s" % CONFIG_NAME)
		
		if os.path.exists(config_file):
			new_config = config_file + ".tmp"
			print "INFO: Configfile %s already exists!\n      Will create %s, which contains the new config.\n      If you want to use the newly generated config you can move %s to %s, otherwise you can delete %s.\n" % (config_file, new_config, new_config, config_file, new_config)

			f =  open(new_config, "wb")
			cfg.write(f)
		else: 
			print "Writing configuration file %s" % config_file
			f =  open(config_file, "wb")
			cfg.write(f)
		
		print "Edit %s to further configure GimmeMotifs." % config_file

module1 = Extension('gimmemotifs.c_metrics', sources = ['gimmemotifs/c_metrics.c'], libraries = ['gsl', 'gslcblas'])


setup (name = 'gimmemotifs',
		cmdclass={"install_data":install_plus_config, "build":custom_build},
		version = '0.50',
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



