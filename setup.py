from setuptools import setup
from setuptools import Extension, Command
from setuptools.command.install import install

from distutils.command.build import build
from distutils.command.install import INSTALL_SCHEMES
from distutils.util import get_platform
from distutils import log as dlog
import distutils.sysconfig
from subprocess import Popen
from platform import machine
from gimmemotifs.tools import *
from gimmemotifs.config import *
from gimmemotifs.shutils import which
from glob import glob
import os
import sys
import shutil
from stat import ST_MODE
import time
import inspect
from io import open

CONFIG_NAME = "gimmemotifs.cfg" 
DESCRIPTION  = """GimmeMotifs is a motif prediction pipeline. 
"""

# trick to get rst file for PyPi, see:
# http://stackoverflow.com/questions/26737222/pypi-description-markdown-doesnt-work/26737672#26737672
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')    
except(IOError, ImportError, RuntimeError):
    long_description = open('README.md', 'r').read()

# are we in the conda build environment?
conda_build = os.environ.get("CONDA_BUILD")

DEFAULT_PARAMS = {
    "max_time": "",
    "analysis": "medium",
    "fraction": "0.2",
    "abs_max": "1000",
    "width": "200",
    "lwidth": "500",
    "pvalue": "0.001",
    "enrichment": "1.5",
    "background": "gc,random",
    "genome": "hg19",
    "tools": "MDmodule,Weeder,MotifSampler",
    "available_tools": "Weeder,MDmodule,MotifSampler,GADEM,MEME,MEMEW,trawler,Improbizer,BioProspector,AMD,ChIPMunk,Jaspar,Homer,XXmotif",
    "cluster_threshold": "0.95",
    "use_strand": "False",
    "markov_model": "1",
    "motif_db": "gimme.vertebrate.v3.1.pwm",
    "scan_cutoff": "0.9",
    "ncpus": "2",
    "use_cache": "False",
}

MOTIF_CLASSES = ["MDmodule", "Meme", "MemeW", "Weeder", "Gadem", "MotifSampler", "Trawler", "Improbizer",  "BioProspector", "Posmo", "ChIPMunk", "Jaspar", "Amd", "Hms", "Homer", "XXmotif"]
LONG_RUNNING = ["GADEM"]


# Included binaries after compile
MOTIF_BINS = {
    "MEME": "src/meme_4.6.0/src/meme.bin",
    "MEMEW": "src/meme_4.6.0/src/meme.bin",
    "MDmodule": "src/MDmodule/MDmodule",
    "BioProspector": "src/BioProspector/BioProspector",
    "Posmo": "src/posmo/posmo",
    "AMD": "src/AMD/AMD.bin",
}

data_files=[
    ('gimmemotifs/templates', 
        [
            'templates/star.png', 
            'templates/report_template.jinja.html', 
            'templates/cluster_template.jinja.html',
            ]),
    
    ('gimmemotifs/templates/sortable', 
        [
            'templates/sortable/sortable.min.js', 
            'templates/sortable/sortable-theme-slick.css',
            ]),
    ('gimmemotifs/score_dists', ['score_dists/total_wic_mean_score_dist.txt']),
    ('gimmemotifs/genes', ['genes/hg18.bed', 'genes/hg19.bed', 'genes/xenTro2.bed', 'genes/mm9.bed']),
    ('gimmemotifs/bg', ['bg/hg19.MotifSampler.bg', 'bg/hg18.MotifSampler.bg', 'bg/mm9.MotifSampler.bg', 'bg/xenTro2.MotifSampler.bg']),
    ('gimmemotifs/motif_databases', [
                                    'motif_databases/JASPAR2010_vertebrate.pwm',
                                    'motif_databases/JASPAR2016_vertebrate.pwm',
                                    'motif_databases/vertebrate_motifs.pwm',
                                    'motif_databases/vertebrate_clusters.pwm',
                                    'motif_databases/vertebrate_clusters.pwm',
                                    'motif_databases/gimme.invertebrate.v1.0.pwm',
                                    'motif_databases/gimme.plant.v1.0.pwm',
                                    'motif_databases/gimme.vertebrate.v3.1.pwm',
                                    'motif_databases/gimme.vertebrate.v3.1.factor2motifs.txt',
                                    'motif_databases/gimme.vertebrate.v3.1.motif2factors.txt',
                                    ]),
#    ('gimmemotifs/doc', ['doc/gimmemotifs_manual.pdf','doc/gimmemotifs_manual.html']),
    ('gimmemotifs/examples', ['examples/TAp73alpha.bed','examples/TAp73alpha.fa']),
    ('gimmemotifs/genome_index', ['genome_index/README.txt'])
]


# Fix for install_data, add share to prefix (borrowed from Dan Christiansen) 
for platform, scheme in INSTALL_SCHEMES.items():
    if platform.startswith('unix_'):
        if scheme['data'][0] == '$' and '/' not in scheme['data']:
            scheme['data'] = os.path.join(scheme['data'], 'share')

class build_tools(Command):
    description = "compile all included motif prediction tools"

    def initialize_options(self):
        self.build_base = None
        self.plat_name = None
        self.build_tools_dir = None
        self.machine = None

    def finalize_options(self):    
        if self.plat_name is None:
            self.plat_name = get_platform()
        self.set_undefined_options('build',('build_base', 'build_base'))
        plat_specifier = ".%s-%s" % (self.plat_name, sys.version[0:3])
        self.build_tools_dir = os.path.join(self.build_base, 'tools' + plat_specifier)
        self.set_undefined_options('build',('custom_build', 'build_tools_dir'))
        self.machine = machine()

    def run(self):
        from compile_externals import compile_all
        prefix = distutils.sysconfig.get_config_var("prefix")
        
        if not os.path.exists(self.build_tools_dir):
            os.mkdir(self.build_tools_dir)

        # Try to compile everything
        compile_all(os.path.join(prefix, "share/gimmemotifs"))

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

        # Copy posmo deps
        if os.path.exists("src/posmo"):
            shutil.copy("src/posmo/clusterwd", self.build_tools_dir)
        
        # Copy ChIPMunk
        if os.path.exists("src/ChIPMunk"):
            if os.path.exists(os.path.join(self.build_tools_dir, "ChIPMunk")):
                shutil.rmtree(os.path.join(self.build_tools_dir, "ChIPMunk"))
            shutil.copytree("src/ChIPMunk", os.path.join(self.build_tools_dir, "ChIPMunk"))
        # Copy HMS
        if os.path.exists("src/HMS"):
            if os.path.exists(os.path.join(self.build_tools_dir, "HMS")):
                shutil.rmtree(os.path.join(self.build_tools_dir, "HMS"))
            shutil.copytree("src/HMS", os.path.join(self.build_tools_dir, "HMS"))

        if self.machine == "x86_64":
            post_fix = "_x86_64"
        else:
            post_fix = "_i386"
        
        # Copy MotifSampler
        if not conda_build and os.path.exists("src/MotifSampler"):
            dlog.info("copying MotifSampler")
            shutil.copy(
                "src/MotifSampler/MotifSampler{}".format(post_fix), 
                    os.path.join(
                        self.build_tools_dir, 
                        "MotifSampler")
                    )
            shutil.copy(
                "src/MotifSampler/CreateBackgroundModel{}".format(post_fix), 
                    os.path.join(
                        self.build_tools_dir, 
                        "CreateBackgroundModel")
                    )
        
        # Copy Improbizer (ameme)
        if os.path.exists("src/Improbizer"):
            dlog.info("copying Improbizer (ameme)")
            shutil.copy("src/Improbizer/ameme{}".format(post_fix), 
                    os.path.join(self.build_tools_dir, "ameme"))

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
        cfg = MotifConfig(use_config="cfg/gimmemotifs.cfg.base")
        
        dlog.info("locating motif programs")
        available = []
        for program in MOTIF_CLASSES:
            # Get class
            m = eval(program)()
            cmd = m.cmd
            
            ### ugly, fixme :)
            if cmd == "ChIPMunk.sh":
                cmd = "ChIPMunk/ChIPMunk.sh"
            if cmd == "hms":
                cmd = "HMS/hms"

            bin = ""
            if cmd == "/bin/false":
                # motif db
                bin = "/bin/false"    
            elif os.path.exists(os.path.join(self.build_tools_dir, cmd)):
                bin = os.path.join(self.build_tools_dir, cmd)
                dlog.info("using included version of %s: %s" % (program, bin))
            else:
                ### ugly, fixme :)
                if     cmd == "ChIPMunk/ChIPMunk.sh":
                    cmd = "ChIPMunk.sh"
                if cmd == "HMS/hms":
                    cmd = "hms"

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
                elif program == "ChIPMunk":
                    dir = bin.replace("ChIPMunk.sh", "")

                available.append(m.name)
                cfg.set_program(m.name, 
                        {
                            "bin":os.path.abspath(bin), 
                            "dir":os.path.abspath(dir),
                            }
                        )

        # Weblogo
        bin = ""
        seq_included = os.path.abspath(os.path.join(self.build_tools_dir, "seqlogo"))
        if os.path.exists(seq_included):
            bin = seq_included
            dlog.info("using included version of weblogo: %s" % seq_included)
        else:
            bin = which("seqlogo")
            dlog.info("using installed version of seqlogo: %s" % (os.path.abspath(bin)))
        if bin:
            cfg.set_seqlogo(os.path.abspath(bin))
        else:
            dlog.info("couldn't find seqlogo")
        
        # Set the available tools in the config file
        DEFAULT_PARAMS["available_tools"] = ",".join(available)
        
        for tool in available:
            if tool in LONG_RUNNING:
                dlog.info("PLEASE NOTE: %s can take a very long time to run on large datasets. Therefore it is not added to the default tools. You can always enable it later, see documentation for details" % tool)
                available.remove(tool)

        DEFAULT_PARAMS["tools"] = ",".join(available)
        cfg.set_default_params(DEFAULT_PARAMS)

        # Write (temporary) config file
        config_file = os.path.join(self.build_cfg, "%s" % CONFIG_NAME)
        dlog.info("writing (temporary) configuration file: %s" % config_file)
        f = open(config_file, "w")
        cfg.write(f)
        f.close()

        # TODO: fix this hack
        my_cfg = open(config_file).read()
        with open(config_file, "w") as f:
            cwd = os.getcwd()
            f.write(my_cfg.replace("/usr/share/gimmemotifs/", cwd + "/"))

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
        for file in self.outfiles:
            #trawler pl's
            if file.endswith("pl"):
                os.chmod(file, 0o755)
    
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
        dlog.info("data_dir: {}".format(data_dir))
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
                bin = cfg.bin(m.name).replace(os.path.abspath(self.build_tools_dir), final_tools_dir) 
                dir = cfg.dir(m.name)
                if dir:
                    dir = dir.replace(os.path.abspath(self.build_tools_dir), final_tools_dir)
                cfg.set_program(m.name, {"bin":bin, "dir":dir})
            
        dir = cfg.get_seqlogo()
        dir = dir.replace(os.path.abspath(self.build_tools_dir), final_tools_dir)
        cfg.set_seqlogo(dir)

        # Use a user-specific configfile if any other installation scheme is used
#        if os.path.abspath(self.install_dir) == "/usr/share":
        config_file = os.path.join(self.install_dir, "gimmemotifs/%s" % CONFIG_NAME)
        self.outfiles = [config_file] 

        if os.path.exists(config_file):
            timestr = time.strftime("%Y%m%d-%H%M%S")        
            old_config = "{}.{}".format(config_file, timestr)
            shutil.move(config_file, old_config)
            dlog.info("INFO: Configfile %s already existed!", config_file)
            dlog.info("INFO: This config has been saved as %s", old_config)
         
        dlog.info("writing configuration file %s" % config_file)
        f =  open(config_file, "w")
        cfg.write(f)
        
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
    
module1 = Extension('gimmemotifs.c_metrics', sources = ['gimmemotifs/c_metrics.c'])

setup (
        name = 'gimmemotifs',
        cmdclass={"build":custom_build, 
                            "build_tools":build_tools,
                            "build_config":build_config,
                            "install":custom_install, 
                            "install_tools":install_tools,
                            "install_config":install_config,
                            },
        version = GM_VERSION,
        long_description = long_description,
        description = DESCRIPTION,
        author = 'Simon van Heeringen',
        author_email = 'simon.vanheeringen@gmail.com',
        url = 'https://github.com/simonvh/gimmemotifs/',
        download_url = 'https://github.com/simonvh/gimmemotifs/tarball/' + GM_VERSION,
        license = 'MIT',
        packages=['gimmemotifs', 'gimmemotifs/commands'],
        ext_modules = [module1],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            ],
        scripts=[
            'scripts/track2fasta.py',
            'scripts/gimme',
            ],
        data_files=data_files,
        install_requires = [
            "setuptools >= 0.7",
            "numpy >= 1.6.0",
            "scipy >= 0.9.0",
            "matplotlib >= 2",
            "jinja2",
            "pyyaml >= 3.10",
            "pybedtools",
            "statsmodels",
            "scikit-learn",
            "sklearn-contrib-lightning",
            "seaborn",
            "pysam",
            "xgboost",
            "xdg",
            "diskcache",
            "xxhash",
            "configparser",
            "six",
            "future",
        ],
)
