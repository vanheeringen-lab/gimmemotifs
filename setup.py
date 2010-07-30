from distutils.core import setup, Extension
from distutils.command.install import install, INSTALL_SCHEMES
from gimmemotifs.utils import which
from gimmemotifs.tools import *
from glob import glob
import os
import sys

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

class install_plus_config(install):
	def run(self):
		install.run(self)
		self.write_config()

	def write_config(self):
		from gimmemotifs.config import MotifConfig
		cfg = MotifConfig(use_config="cfg/gimmemotifs.cfg.example")

		data_dir = os.path.abspath(self.install_data)
		cfg.set_template_dir(os.path.join(data_dir, 'gimmemotifs/templates'))
		cfg.set_gene_dir(os.path.join(data_dir, 'gimmemotifs/genes'))
		cfg.set_score_dir(os.path.join(data_dir, 'gimmemotifs/score_dists'))
		cfg.set_index_dir(os.path.join(data_dir, 'gimmemotifs/genome_index'))
		cfg.set_motif_dir(os.path.join(data_dir, 'gimmemotifs/motif_databases'))
		cfg.set_bg_dir(os.path.join(data_dir, 'gimmemotifs/bg'))
		
		print 
		print "Trying to locate motif programs"
		
		MOTIF_CLASSES = ["MDmodule", "Meme", "Weeder", "Gadem", "MotifSampler", "Trawler", "Improbizer", "MoAn", "BioProspector"]
		available = []
		for program in MOTIF_CLASSES:
			m = eval(program)()
			cmd = m.cmd
			bin = which(cmd)
			if bin:
				print "Found %s in %s" % (m.name, bin)
				available.append(m.name)
				dir = None
				if program == "Weeder":
					dir = bin.replace("weederTFBS.out","")
				elif program == "Meme":
					dir = bin.replace("bin/meme", "")
				elif program == "Trawler":
					dir = bin.replace("bin/trawler.pl", "")
		
				cfg.set_program(m.name, {"bin":bin, "dir":dir})
			else:
				print "Couldn't find %s" % m.name
		
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
		if os.path.abspath(self.install_data) == "/usr/share":
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
		cmdclass={"install":install_plus_config},
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
		data_files=[
			('gimmemotifs/templates', ['templates/cluster_template.kid', 'templates/report_template.kid', 'templates/report_template_v2.kid', 'templates/cluster_template_v2.kid']),
			('gimmemotifs/score_dists', ['score_dists/total_wic_mean_score_dist.txt']),
			('gimmemotifs/genes', ['genes/hg18.bed', 'genes/mm9.bed', 'genes/xenTro2.bed']),
			('gimmemotifs/bg', ['bg/hg18.MotifSampler.bg', 'bg/mm9.MotifSampler.bg', 'bg/xenTro2.MotifSampler.bg']),
			('gimmemotifs/motif_databases', ['motif_databases/jaspar.pfm']),
			('gimmemotifs/genome_index', ['genome_index/README.txt'])
			]
)



