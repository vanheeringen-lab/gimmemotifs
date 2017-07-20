# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""Interface module for all motif programs."""
from __future__ import print_function
# Python imports
import re
import os
import sys
from subprocess import Popen, PIPE
import shutil
from tempfile import NamedTemporaryFile, mkdtemp
import io
import glob

# gimme imports
from gimmemotifs.config import MotifConfig
from gimmemotifs.shutils import which

# Necessary for when setup.py needs to import 
# refactor at some point
try:
    from gimmemotifs.motif import read_motifs, Motif
    from gimmemotifs.fasta import Fasta
except ImportError:
    pass

def get_tool(name): 
    """
    Returns an instance of a specific tool.

    Parameters
    ----------
    name : str
        Name of the tool (case-insensitive).

    Returns
    -------
    tools : MotifProgram instance
    """
    tool = name.lower()
    if tool not in __tools__:
        raise ValueError("Tool {0} not found!\n".format(name))

    t = __tools__[tool]()

    if not t.is_installed():
        sys.stderr.write("Tool {0} not installed!\n".format(tool))

    if not t.is_configured():
        sys.stderr.write("Tool {0} not configured!\n".format(tool))

    return t

def locate_tool(name, verbose=True): 
    """
    Returns the binary of a tool.

    Parameters
    ----------
    name : str
        Name of the tool (case-insensitive).

    Returns
    -------
    tool_bin : str
        Binary of tool.
    """
    m = get_tool(name) 
    tool_bin = which(m.cmd) 
    if tool_bin:
        if verbose:
            print("Found {} in {}".format(m.name, tool_bin)) 
        return tool_bin 
    else: 
        print("Couldn't find {}".format(m.name))

class MotifProgram(object):
    
    """Motif program base class."""

    config = MotifConfig()
    local_bin = None

    def __init__(self):
        pass

    def bin(self):
        """
        Get the command used to run the tool.

        Returns
        -------
        command : str
            The tool system command.
        """
        if self.local_bin:
            return self.local_bin
        else:
            return self.config.bin(self.name)

    def dir(self):
        """
        Get the installation directory of the tool.

        Returns
        -------
        dir : str
            The tool directory.
        """
        return self.config.dir(self.name)

    def is_configured(self):
        """
        Check if the tool is configured.
        
        Returns
        -------
        is_configured : bool
            True if the tool is configured.
        """
        return self.config.is_configured(self.name)
    
    def is_installed(self):
        """
        Check if the tool is installed.
        
        Returns
        -------
        is_installed : bool
            True if the tool is installed.
        """
        return self.is_configured() and os.access(self.bin(), os.X_OK)

    def run(self, fastafile, params=None, tmp=None):
        """
        Run the tool and predict motifs from a FASTA file.

        Parameters
        ----------
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        tmp : str, optional
            Directory to use for creation of temporary files.
       
        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        if not self.is_configured():
            raise ValueError("%s is not configured" % self.name)

        if not self.is_installed():
            raise ValueError("%s is not installed or not correctly configured" % self.name)
        
        self.tmpdir = mkdtemp(prefix="{0}.".format(self.name), dir=tmp)
        fastafile = os.path.abspath(fastafile)
 
        try:
            return self._run_program(self.bin(), fastafile, params)
        except KeyboardInterrupt:
            return ([], "Killed", "Killed")
#        except Exception as e:
#            return ([], "", e.strerror)

class XXmotif(MotifProgram):

    """
    Predict motifs using XXmotif.
    
    Reference: 
    """

    def __init__(self):
        self.name = "XXmotif"
        self.cmd = "XXmotif"
        self.use_width = False
        self.default_params = {
                "single":False, 
                "background":None, 
                "analysis":"medium", 
                "number":5, 
                "width":10,
                }

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        if prm["background"]:
            # Absolute path, just to be sure
            prm["background"] =  os.path.abspath(prm["background"])
            prm["background"] = " --negSet {0} ".format(
                    prm["background"])
        
        prm["strand"] = ""
        if not prm["single"]:
            prm["strand"] = " --revcomp "

        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run XXmotif and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params) 
        
        outfile = os.path.join(
                self.tmpdir, 
                os.path.basename(fastafile.replace(".fa", ".pwm")))
        
        stdout = ""
        stderr = ""
        
        cmd = "%s %s %s --localization --batch %s %s" % (
            bin,
            self.tmpdir, 
            fastafile,
            params["background"],
            params["strand"],
            )
        
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()
        
        motifs = []
        
        if os.path.exists(outfile):
            with open(outfile) as f:
                motifs = read_motifs(f, fmt="xxmotif")
            for m in motifs:
                m.id = "{0}_{1}".format(self.name, m.id)
        else:
            stdout += "\nMotif file {0} not found!\n".format(outfile)
            stderr += "\nMotif file {0} not found!\n".format(outfile)
        
        return motifs, stdout, stderr
        
class Homer(MotifProgram):

    """
    Predict motifs using Homer.
    
    Reference: Heinz et al, 2010; PMID: 20513432 
    """

    def __init__(self):
        self.name = "Homer"
        self.cmd = "homer2"
        self.use_width = True
        self.default_params = {
                "single":False, 
                "background":None, 
                "analysis":"medium", 
                "number":5, "width":10
                }

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Background file is essential!
        if not prm["background"]:
            print("Background file needed!")
            sys.exit()
        
        prm["background"] =  os.path.abspath(prm["background"])
        
        prm["strand"] = ""
        if prm["single"]:
            prm["strand"] = " -strand + "
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Homer and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params) 
        
        outfile = NamedTemporaryFile(
                mode="w",
                dir=self.tmpdir, 
                prefix= "homer_w{}.".format(params["width"])
                ).name
        
        cmd = "%s denovo -i %s -b %s -len %s -S %s %s -o %s -p 8" % (
            bin,
            fastafile,
            params["background"],
            params["width"],
            params["number"],
            params["strand"],
            outfile)

        stderr = ""
        stdout = "Running command:\n{}\n".format(cmd)
        
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()
        
        motifs = []
        
        if os.path.exists(outfile):
            with open(outfile) as f:
                motifs = read_motifs(f, fmt="pwm")
                for i, m in enumerate(motifs):
                    m.id = "{}_{}_{}".format(self.name, params["width"], i + 1)
        
        return motifs, stdout, stderr

class BioProspector(MotifProgram):

    """
    Predict motifs using BioProspector.
    
    Reference: 
    """

    def __init__(self):
        self.name = "BioProspector"
        self.cmd = "BioProspector"
        self.use_width = True
        self.default_params = {
                "single":False, 
                "background":None, 
                "analysis":"medium", 
                "number":5, 
                "width":10}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
        
        # Background file is essential!
        if not prm["background"]:
            print("Background file needed!")
            sys.exit()
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        prm["strand"] = 2
        if prm["single"]:
            prm["strand"] = 1
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run BioProspector and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)
        
        outfile = os.path.join(self.tmpdir, "bioprospector.out")    
        
        stdout = ""
        stderr = ""
        
        cmd = "%s -i %s -W %s -d %s -b %s -r %s -o %s" % (
            bin,
            fastafile,
            params["width"],
            params["strand"],
            params["background"],
            params["number"],
            outfile)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()
        
        motifs = []
        
        if os.path.exists(outfile):
            with open(outfile) as f: 
                motifs = self.parse(f)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert BioProspector output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing BioProspector output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        
        p = re.compile(r'^\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
        pwm = []
        motif_id = ""
        for line in fo.readlines():
            if line.startswith("Motif #"):
                if pwm:
                    m = Motif(pwm)
                    m.id = "BioProspector_w%s_%s" % (len(m), motif_id)
                    motifs.append(m)
                motif_id =  line.split("#")[1].split(":")[0]
                pwm = []
            else:
                m = p.search(line)
                if m:
                    pwm.append([float(m.group(x))/100.0 for x in range(1,5)])

        if pwm:
            m = Motif(pwm)
            m.id = "BioProspector_w%s_%s" % (len(m), motif_id)
            motifs.append(m)
        return motifs


class Hms(MotifProgram):
    
    """
    Predict motifs using HMS.
    
    Reference: 
    """

    def __init__(self):
        self.name = "HMS"
        self.cmd = "hms"
        self.use_width = True
        self.default_params = {"background":None}
     
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 
  
    def _prepare_files(self, fastafile):

        hmsdir = self.dir() 
        thetas = ["theta%s.txt" % i for i in [0,1,2,3]]
        for t in thetas:
            shutil.copy(os.path.join(hmsdir, t), self.tmpdir)

        summitfile = os.path.join(self.tmpdir, "HMS.in.summits.txt")
        outfile = os.path.join(self.tmpdir, "thetafinal.txt")    
        fgfile = os.path.join(self.tmpdir, "HMS.in.fa")
        
        shutil.copy(fastafile, fgfile)
        fa = Fasta(fgfile)
        with open(summitfile, "w") as out:
            for seq in fa.seqs:
                out.write("%s\n" % (len(seq) / 2))
        return fgfile, summitfile, outfile 

    
    def _run_program(self, bin, fastafile, params=None):
        """
        Run HMS and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)
        
        default_params = {"width":10}
        if params is not None: 
            default_params.update(params)
        
        fgfile, summitfile, outfile = self._prepare_files(fastafile)
                
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        cmd = "{} -i {} -w {} -dna 4 -iteration 50 -chain 20 -seqprop -0.1 -strand 2 -peaklocation {} -t_dof 3 -dep 2".format(
                bin, 
                fgfile, 
                params['width'], 
                summitfile)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        stdout,stderr = p.communicate()
        
        os.chdir(current_path)
        
        motifs = []
        if os.path.exists(outfile):
            with open(outfile) as f: 
                motifs = self.parse(f)
                for i,m in enumerate(motifs):
                    m.id = "HMS_w{}_{}".format(params['width'], i + 1)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert HMS output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing HMS output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        m = [[float(x) for x in fo.readline().strip().split(" ")] for i in range(4)]
        matrix = [[m[0][i], m[1][i],m[2][i],m[3][i]] for i in range(len(m[0]))]
        motifs = [Motif(matrix)]
        motifs[-1].id = self.name
        
        return motifs

class Amd(MotifProgram):
    
    """
    Predict motifs using AMD.

    Reference: 
    """

    def __init__(self):
        self.name = "AMD"
        self.cmd = "AMD.bin"
        self.use_width = False
        self.default_params = {"background":None}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Background file is essential!
        if not prm["background"]:
            raise ValueError("Background file needed!")
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])

        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run AMD and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)

        fgfile = os.path.join(self.tmpdir, "AMD.in.fa")
        outfile = fgfile + ".Matrix"    
        shutil.copy(fastafile, fgfile)
        
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        stdout = ""
        stderr = ""
    
        cmd = "%s -F %s -B %s" % (
                bin, 
                fgfile, 
                params["background"],
                )
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()
        
        os.chdir(current_path)
        motifs = []
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert AMD output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing AMD output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        
        #160:  112  CACGTGC      7.25   chr14:32308489-32308689
        p = re.compile(r'\d+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)')
        wm = []
        name = ""
        for line in fo.readlines():
            if line.startswith("Motif") and line.strip().endswith(":"):
                if name:
                    motifs.append(Motif(wm))
                    motifs[-1].id = name
                    name = ""
                    wm = []
                name = "%s_%s" % (self.name, line.split(":")[0])
            else:
                m = p.search(line)
                if m:
                    wm.append([float(m.group(x)) for x in range(1,5)])
        motifs.append(Motif(wm))
        motifs[-1].id = name
        
        return motifs

class Improbizer(MotifProgram):
    
    """
    Predict motifs using Improbizer.

    Reference: 
    """

    def __init__(self):
        self.name = "Improbizer"
        self.cmd = "ameme"
        self.use_width = False 
        self.default_params = {"background":None, "number":10}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
        
        # Not strictly necessary, but recommended
        if not params["background"]:
            print("Background file needed!")
            sys.exit()
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        prm["outfile"] = os.path.join(self.tmpdir, "improbizer.out.html")    
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Improbizer and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params) 
        
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        stdout = ""
        stderr = ""
        cmd = "%s good=%s bad=%s numMotifs=%s > %s" % (
                bin, 
                fastafile, 
                params["background"],
                params["number"],
                params["outfile"])
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()
        
        motifs = []
        if os.path.exists(params["outfile"]):
            f = open(params["outfile"])
            motifs = self.parse(f)
            f.close()
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert Improbizer output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing Improbizer output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        p = re.compile(r'\d+\s+@\s+\d+\.\d+\s+sd\s+\d+\.\d+\s+(\w+)$')

        line = fo.readline()
        while line and line.find("Color") == -1:
            m = p.search(line)
            if m:
                pwm_data = {}
                for i in range(4):
                    vals = [x.strip() for x in fo.readline().strip().split(" ") if x]
                    pwm_data[vals[0].upper()] = vals[1:]
                pwm = []
                for i in range(len(pwm_data["A"])):
                    pwm.append([float(pwm_data[x][i]) for x in ["A","C","G","T"]])
                motifs.append(Motif(pwm))
                motifs[-1].id = "%s_%s" % (self.name, m.group(1))
            line = fo.readline()
        
        return motifs

class Trawler(MotifProgram):
    
    """
    Predict motifs using Trawler.

    Reference: Ettwiller, 2010; PMID: 17589518 
    """

    def __init__(self):
        self.name = "trawler"
        self.cmd = "trawler"
        self.use_width = False
        self.default_params = {"single":False, "background":None}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Background file is essential!
        if not prm["background"]:
            print("Background file needed!")
            sys.exit()
        
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
         
        prm['strand'] = "double"
        if prm["single"]:
            prm["strand"] = "single"
       
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Trawler and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)

        tmp = NamedTemporaryFile(mode="w", dir=self.tmpdir, delete=False)
        shutil.copy(fastafile, tmp.name)
        fastafile = tmp.name
    
        current_path = os.getcwd()
        os.chdir(self.dir())
        
        motifs = []
        stdout = ""
        stderr = ""
        for wildcard in [0,1,2]:
            cmd = "%s -sample %s -background %s -directory %s -strand %s -wildcard %s" % (
                    bin, 
                    fastafile, 
                    params["background"], 
                    self.tmpdir, 
                    params["strand"],
                    wildcard,
                    )
            
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
            out,err = p.communicate()
            stdout += out.decode()
            stderr += err.decode()
    
            os.chdir(current_path)
            pwmfiles = glob.glob("{}/tmp*/result/*pwm".format(self.tmpdir))
            if len(pwmfiles) > 0:
                out_file = pwmfiles[0]
                stdout += "\nOutfile: {}".format(out_file)
                 
                my_motifs = []
                if os.path.exists(out_file):
                    with open(out_file) as f: 
                        my_motifs = read_motifs(f, fmt="pwm")
                    for m in motifs:
                        m.id = "{}_{}".format(self.name, m.id)
                    stdout += "\nTrawler: {} motifs".format(len(motifs))
            
                # remove temporary files
                if os.path.exists(tmp.name):
                    os.unlink(tmp.name)
            
                for motif in my_motifs:
                    motif.id = "{}_{}_{}".format(self.name, wildcard, motif.id)
            
                motifs += my_motifs
            else:
                stderr += "\nNo outfile found"

        return motifs, stdout, stderr

class Weeder(MotifProgram):
    
    """
    Predict motifs using Weeder.
    
    Reference: 
    
    """
    
    def __init__(self):
        self.name = "Weeder"
        self.cmd = "weeder2"
        self.use_width = False
        self.default_params = {"organism":"hg19", "single":False}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 
    
    def _run_program(self, bin,fastafile, params=None):
        """
        Run Weeder and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)
        
        organism = params["organism"]
        weeder_organisms = {
            "hg18":"HS",
            "hg19":"HS",
            "hg38":"HS",
            "mm9":"MM",
            "mm10":"MM",
            "dm3":"DM",
            "dm5":"DM",
            "dm6":"DM",
            "yeast":"SC",
            "sacCer2":"SC",
            "sacCer3":"SC",
            "TAIR10":"AT",
            "TAIR11":"AT",
            }
        weeder_organism = weeder_organisms.get(organism, "HS")

        tmp = NamedTemporaryFile(dir=self.tmpdir)
        name = tmp.name
        tmp.close()
        shutil.copy(fastafile, name)
        fastafile = name
    
        cmd = "{} -f {} -O".format(
                self.cmd, 
                fastafile,
                weeder_organism,
                )
        
        if params["single"]:
            cmd += " -ss"
        
        #print cmd
        stdout, stderr = "", ""
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir) 
        out,err = p.communicate()
        stdout += out.decode()
        stderr += err.decode()

        motifs = []
        if os.path.exists(fastafile + ".matrix.w2"):
            f = open(fastafile + ".matrix.w2")
            motifs = self.parse(f)
            f.close()
        
        for m in motifs:
            m.id = "{}_{}".format(self.name, m.id.split("\t")[0])
        
        for ext in [".w2", ".matrix.w2" ]:
            if os.path.exists(fastafile + ext):
                os.unlink(fastafile + ext)

        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert Weeder output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing Weeder output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        return read_motifs(fo, fmt="jaspar")

class MotifSampler(MotifProgram):
    
    """
    Predict motifs using MotifSampler.
    
    Reference: 
    """
    
    def __init__(self):
        self.name = "MotifSampler"
        self.cmd = "MotifSampler"
        self.use_width = True
        self.default_params = {
                "width":10, 
                "background_model":"", 
                "single":False, 
                "number":10}
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        if prm["background_model"]:
            # Absolute path, just to be sure
            prm["background_model"] = os.path.abspath(prm["background_model"])
        else:
            if prm.get("organism", None):
                prm["background_model"] = os.path.join(
                        self.config.get_bg_dir(), 
                        "{}.{}.bg".format(
                            prm["organism"], 
                            "MotifSampler"))
            else:            
                raise Exception("No background specified for {}".format(self.name))
        
        prm["strand"] = 1
        if prm["single"]:
            prm["strand"] = 0
        
        tmp = NamedTemporaryFile(dir=self.tmpdir)
        prm["pwmfile"] = tmp.name

        tmp2  = NamedTemporaryFile(dir=self.tmpdir)
        prm["outfile"] = tmp2.name
 
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run MotifSampler and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)
        # TODO: test organism
        #cmd = "%s -f %s -b %s -m %s -w %s -n %s -o %s -s %s > /dev/null 2>&1" % (
        cmd = "%s -f %s -b %s -m %s -w %s -n %s -o %s -s %s" % (
                bin, 
                fastafile, 
                params["background_model"], 
                params["pwmfile"], 
                params["width"], 
                params["number"], 
                params["outfile"],
                params["strand"],
                )
        #print cmd
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        stdout, stderr = p.communicate()
        
        #stdout,stderr = "",""
        #p = Popen(cmd, shell=True)
        #p.wait()

        motifs = []
        if os.path.exists(params["outfile"]):
            with open(params["outfile"]) as f:
                motifs = self.parse_out(f)
        
        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MotifSampler output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing MotifSampler output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        pwm = []
        info = {}
        for line in fo.readlines():
            if line.startswith("#"):
                vals =  line.strip()[1:].split(" = ")
                if len(vals) > 1:
                    info[vals[0]] = vals[1]
            elif len(line) > 1:
                pwm.append([float(x) for x in line.strip().split("\t")])
            else:
                motifs.append(Motif())
                motifs[-1].consensus = info["Consensus"]
                motifs[-1].width = info["W"]
                motifs[-1].id = info["ID"]
                motifs[-1].pwm = pwm[:]
                pwm = []
            
        return motifs

    def parse_out(self, fo):
        """
        Convert MotifSampler output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing MotifSampler output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}
        pseudo = 0.0 # Should be 1/sqrt(# of seqs)
        aligns = {}
        for line in fo.readlines():
            if line.startswith("#"):
                pass
            elif len(line) > 1:
                vals = line.strip().split("\t")
                m_id, site = [x.strip().split(" ")[1].replace('"',"") for x in vals[8].split(";") if x]
                #if vals[6] == "+":
                if site.upper().find("N") == -1:
                    aligns.setdefault(m_id, []).append(site)
                #else:
                #    print site, rc(site)
                #    aligns.setdefault(id, []).append(rc(site))
                        
        for m_id, align in aligns.items():
            #print id, len(align)

            width = len(align[0])
            pfm =  [[0 for x in range(4)] for x in range(width)]
            for row in align:
                for i in range(len(row)):
                    pfm[i][nucs[row[i]]] += 1
            total = float(len(align))
            pwm = [[(x + pseudo/4)/total+(pseudo) for x in row] for row in pfm]
            m = Motif()
            m.align = align[:]
            m.pwm = pwm[:]
            m.pfm = pfm[:]
            m.id = m_id
            motifs.append(m)    
        return motifs

class MDmodule(MotifProgram):
    
    """
    Predict motifs using MDmodule.

    Reference: 
    """

    def __init__(self):
        self.name = "MDmodule"
        self.cmd = "MDmodule"
        self.use_width = True
         
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 
   
    def _run_program(self, bin, fastafile, params=None):
        """
        Run MDmodule and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        default_params = {"width":10, "number":10}
        if params is not None: 
            default_params.update(params)
        
        new_file = os.path.join(self.tmpdir, "mdmodule_in.fa")
        shutil.copy(fastafile, new_file)
        
        fastafile = new_file
        pwmfile = fastafile + ".out"
    
        width = default_params['width']
        number = default_params['number']
    
        current_path = os.getcwd()
        os.chdir(self.tmpdir)    
        cmd = "%s -i %s -a 1 -o %s -w %s -t 100 -r %s" % (bin, fastafile, pwmfile, width, number)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        stdout,stderr = p.communicate()
        
        stdout = "cmd: {}\n".format(cmd) + stdout.decode() 
            
        motifs = []
        if os.path.exists(pwmfile):
            with open(pwmfile) as f:
                motifs = self.parse(f)
        
        os.chdir(current_path)
        
        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MDmodule output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing MDmodule output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}
        p = re.compile(r'(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
        pf = re.compile(r'>.+\s+[bf]\d+\s+(\w+)')

        pwm = []
        pfm = []
        align = []
        m_id = ""
        for line in fo.readlines():
            if line.startswith("Motif"):
                if m_id:
                    motifs.append(Motif())
                    motifs[-1].id = m_id
                    motifs[-1].pwm = pwm
                    motifs[-1].pfm = pfm
                    motifs[-1].align = align
                    pwm = []
                    pfm = []
                    align = []
                m_id = line.split("\t")[0]
            else: 
                m = p.search(line)
                if m:
                    pwm.append([float(m.group(x))/100 for x in [2,3,4,5]])
                m = pf.search(line)
                if m:
                    if not pfm:
                        pfm = [[0 for x in range(4)] for x in range(len(m.group(1)))]
                    for i in range(len(m.group(1))):
                        pfm[i][nucs[m.group(1)[i]]] += 1
                            
                    align.append(m.group(1))
        
        if pwm:
            motifs.append(Motif())
            motifs[-1].id = m_id
            motifs[-1].pwm = pwm
            motifs[-1].pfm = pfm
            motifs[-1].align = align

        return motifs

class ChIPMunk(MotifProgram):
    
    """
    Predict motifs using ChIPMunk.

    Reference: 
    """
    
    def __init__(self):
        self.name = "ChIPMunk"
        self.cmd = "ChIPMunk.sh"
        self.use_width = True
        self.default_params = {}

    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run ChIPMunk and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        params = self._parse_params(params)
        basename = "munk_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        out = open(new_file, "w")
        f = Fasta(fastafile)
        for seq in f.seqs:
            header = " ".join(["%0.1f" % x for x in list(range(len(seq) // 2)) + list(range(len(seq) // 2, 0, -1))])
            out.write(">%s\n" % header)
            out.write("%s\n" % seq)
        out.close()
        
        fastafile = new_file
        outfile = fastafile + ".out"

        current_path = os.getcwd()
        os.chdir(self.dir())
       
        # Max recommended by ChIPMunk userguide
        ncpus = 4
        cmd = "{} {} {} y, 1.0 s:{} 100 10 1 {} 1>{}".format(
                bin, 
                params.get("width", 8),
                params.get("width", 20),
                fastafile, 
                ncpus, 
                outfile
                )
        p = Popen(cmd, shell=True, stderr=PIPE) 
        stdout, stderr = p.communicate()
        if "RuntimeException" in stderr.decode():
            return [], stdout, stderr
        
        motifs = []
        if os.path.exists(outfile):
            with open(outfile) as f:
                motifs = self.parse(f)
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr
        
    def parse(self, fo):
        """
        Convert ChIPMunk output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing ChIPMunk output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        #KDIC|6.124756232026243
        #A|517.9999999999999 42.99999999999999 345.99999999999994 25.999999999999996 602.9999999999999 155.99999999999997 2.9999999999999996 91.99999999999999
        #C|5.999999999999999 4.999999999999999 2.9999999999999996 956.9999999999999 91.99999999999999 17.999999999999996 22.999999999999996 275.99999999999994
        #G|340.99999999999994 943.9999999999999 630.9999999999999 6.999999999999999 16.999999999999996 48.99999999999999 960.9999999999999 14.999999999999998
        #T|134.99999999999997 7.999999999999999 19.999999999999996 9.999999999999998 287.99999999999994 776.9999999999999 12.999999999999998 616.9999999999999
        #N|999.9999999999998
        line = fo.readline()
        if not line:
            return []
        
        while not line.startswith("A|"):
            line = fo.readline() 
        matrix = []
        for _ in range(4):
            matrix.append([float(x) for x in line.strip().split("|")[1].split(" ")])
            line = fo.readline()
        #print matrix
        matrix = [[matrix[x][y] for x in range(4)] for y in range(len(matrix[0]))]
        #print matrix
        m = Motif(matrix)
        m.id = "ChIPMunk_w%s" % len(m)
        return [m]


class Posmo(MotifProgram):
    
    """
    Predict motifs using Posmo.

    Reference: 
    """
    
    def __init__(self):
        self.name = "Posmo"
        self.cmd = "posmo"
        self.use_width = False
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run Posmo and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        default_params = {}
        if params is not None: 
            default_params.update(params)
        
        basename = "posmo_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        shutil.copy(fastafile, new_file)
        
        fastafile = new_file
        #pwmfile = fastafile + ".pwm"
    
        motifs = []
        current_path = os.getcwd()
        os.chdir(self.tmpdir)    
        for x in ["111111", "11111111"]:
            outfile = "%s.%s.out" % (fastafile, x)
            cmd = "%s 5000 %s %s 1.6 2.5 20 200" % (bin, x, fastafile)
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
            stdout, stderr = p.communicate()
            stdout = stdout.decode()
            stderr = stderr.decode()

            context_file = fastafile.replace(basename, "context.%s.%s.txt" % (basename, x))
            cmd = "%s %s %s simi.txt 0.88 10 2 10" % (bin.replace("posmo","clusterwd"), context_file, outfile)
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
            out, err = p.communicate()
            stdout += out.decode()
            stderr += err.decode()
        
            if os.path.exists(outfile):
                with open(outfile) as f:
                    motifs += self.parse(f)
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert Posmo output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing Posmo output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []

        lines = [fo.readline() for x in range(6)]
        while lines[0]:
            matrix = [[float(x) for x in line.strip().split("\t")] for line in lines[2:]]
            matrix = [[matrix[x][y] for x in range(4)] for y in range(len(matrix[0]))]
            m = Motif(matrix)
            m.id = lines[0].strip().split(" ")[-1]
            motifs.append(m)
            lines = [fo.readline() for x in range(6)]
        
        for i,motif in enumerate(motifs):
            motif.id = "%s_%s" % (self.name, i + 1)
            motif.trim(0.25)
        
        return motifs

class Gadem(MotifProgram):
    
    """
    Predict motifs using GADEM.

    Reference: 
    """
    
    def __init__(self):
        self.name = "GADEM"
        self.cmd = "gadem"
        self.use_width = False
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run GADEM and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        default_params = {}
        if params is not None: 
            default_params.update(params)
        
        new_file = os.path.join(self.tmpdir, "gadem_in.fa")
        shutil.copy(fastafile, new_file)
        
        fastafile = new_file
        pwmfile = fastafile + ".pwm"
        outfile = fastafile + ".out"
    
        current_path = os.getcwd()
        os.chdir(self.tmpdir)    
        cmd = "%s -fseq %s -fpwm %s -fout %s" % (bin, fastafile, pwmfile, outfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        stdout, stderr = p.communicate()
            
        motifs = []
        if os.path.exists(pwmfile):
            with open(pwmfile) as f:
                motifs = self.parse(f)
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert GADEM output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing GADEM output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

        lines = fo.readlines()
        for i in range(0, len(lines), 5):
            align = []
            pwm = []
            pfm = []
            m_id = ""
            line = lines[i].strip()
            m_id = line[1:]
            number = m_id.split("_")[0][1:]
            if os.path.exists("%s.seq" % number):
                with open("%s.seq" % number) as f:
                    for l in f:
                        if "x" not in l and "n" not in l:
                            l = l.strip().upper()
                            align.append(l)
                            if not pfm:
                                pfm = [[0 for x in range(4)] for x in range(len(l))]
                            for p in range(len(l)):
                                pfm[p][nucs[l[p]]] += 1
    
            m = [l.strip().split(" ")[1].split("\t") for l in lines[i + 1: i + 5]]

            pwm = [[float(m[x][y]) for x in range(4)] for y in range(len(m[0]))]


            motifs.append(Motif(pwm))
            motifs[-1].id = "{}_{}".format(self.name, m_id)
            #motifs[-1].pwm = pwm
            if align:
                motifs[-1].pfm = pfm
                motifs[-1].align = align

        return motifs

class Jaspar(MotifProgram):
    def __init__(self):
        self.name = "JASPAR"
        self.cmd = "/bin/false"    
        self.use_width = False
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Get enriched JASPAR motifs in a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        fname = os.path.join(self.config.get_motif_dir(), "JASPAR2010_vertebrate.pwm")
        with open(fname) as f:
            motifs =  read_motifs(f, fmt="pwm")

        for motif in motifs:
            motif.id = "JASPAR_%s" % motif.id
        return motifs, "", ""

class Meme(MotifProgram):
    
    """
    Predict motifs using MEME.
    
    Reference: 
    """

    def __init__(self):
        self.name = "MEME"
        self.cmd = "meme.bin"
        self.use_width = True
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run MEME and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        default_params = {"width":10, "single":False, "number":10}
        if params is not None: 
            default_params.update(params)
        
        tmp = NamedTemporaryFile(dir=self.tmpdir)
        tmpname = tmp.name
    
        strand = "-revcomp"
        width = default_params["width"]
        number = default_params["number"]
        
        cmd = [bin, fastafile, "-text","-dna","-nostatus","-mod", "zoops","-nmotifs", "%s" % number, "-w","%s" % width, "-maxsize", "10000000"]
        if not default_params["single"]:
            cmd.append(strand)
        
        #sys.stderr.write(" ".join(cmd) + "\n")
        p = Popen(cmd, bufsize=1, stderr=PIPE, stdout=PIPE) 
        stdout,stderr = p.communicate()

        motifs = []
        motifs = self.parse(io.StringIO(stdout.decode()))
        
        # Delete temporary files
        tmp.close()
         
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MEME output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing MEME output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

        p = re.compile('BL   MOTIF (\d+) width=(\d+) seqs=(\d+)')
        pa = re.compile('\)\s+(\w+)')
        line = fo.readline()
        while line:
            m = p.search(line)
            align = []
            pfm = []    
            if m:
                id = "%s_%s_w%s" % (self.name, m.group(1), m.group(2))
                while not line.startswith("//"):
                    ma = pa.search(line)
                    if ma:
                        l = ma.group(1)
                        align.append(l)
                        if not pfm:
                            pfm = [[0 for x in range(4)] for x in range(len(l))]
                        for pos in range(len(l)):
                            if l[pos] in nucs:
                                pfm[pos][nucs[l[pos]]] += 1
                            else:
                                for i in range(4):
                                    pfm[pos][i] += 0.25
                    
                    line = fo.readline()
                
                motifs.append(Motif(pfm[:]))
                motifs[-1].id = id
                motifs[-1].align = align[:]
            line = fo.readline()

        return motifs

class MemeW(MotifProgram):
    
    """
    Predict motifs using MEME
    
    Reference: 
    """
    
    def __init__(self):
        self.name = "MEMEW"
        self.cmd = "meme.bin"
        self.use_width = False
    
    def _parse_params(self, params=None):
        """
        Parse parameters.

        Combine default and user-defined parameters.
        """
        prm = self.default_params.copy()
        if params is not None: 
            prm.update(params)
 
        # Absolute path, just to be sure
        prm["background"] =  os.path.abspath(prm["background"])
        
        return prm 

    def _run_program(self, bin, fastafile, params=None):
        """
        Run MEME and predict motifs from a FASTA file.

        Parameters
        ----------
        bin : str
            Command used to run the tool.
        
        fastafile : str
            Name of the FASTA input file.

        params : dict, optional
            Optional parameters. For some of the tools required parameters
            are passed using this dictionary.

        Returns
        -------
        motifs : list of Motif instances
            The predicted motifs.

        stdout : str
            Standard out of the tool.
        
        stderr : str
            Standard error of the tool.
        """
        default_params = {"single":False, "number":5}
        if params is not None: 
            default_params.update(params)
        
        tmp = NamedTemporaryFile(dir=self.tmpdir)
        tmpname = tmp.name
    
        strand = "-revcomp"
        number = default_params["number"]
        
        cmd = [bin, fastafile, "-text","-dna","-nostatus","-mod", "zoops","-nmotifs", "%s" % number, "-minw", "6", "-maxw","20", "-maxsize", "10000000"]
        if not default_params["single"]:
            cmd.append(strand)
        
        #sys.stderr.write(" ".join(cmd) + "\n")
        p = Popen(cmd, bufsize=1, stderr=PIPE, stdout=PIPE) 
        stdout,stderr = p.communicate()

        motifs = []
        motifs = self.parse(io.StringIO(stdout.decode()))
        
        # Delete temporary files
        tmp.close()
         
        return motifs, stdout, stderr

    def parse(self, fo):
        """
        Convert MEME output to motifs
        
        Parameters
        ----------
        fo : file-like
            File object containing MEME output.

        Returns
        -------
        motifs : list
            List of Motif instances.
        """
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

        p = re.compile('BL   MOTIF (\d+) width=(\d+) seqs=(\d+)')
        pa = re.compile('\)\s+(\w+)')
        line = fo.readline()
        while line:
            m = p.search(line)
            align = []
            pfm = []    
            if m:
                id = "%s_%s_w%s" % (self.name, m.group(1), m.group(2))
                while not line.startswith("//"):
                    ma = pa.search(line)
                    if ma:
                        l = ma.group(1)
                        align.append(l)
                        if not pfm:
                            pfm = [[0 for x in range(4)] for x in range(len(l))]
                        for pos in range(len(l)):
                            if l[pos] in nucs:
                                pfm[pos][nucs[l[pos]]] += 1
                            else:
                                for i in range(4):
                                    pfm[pos][i] += 0.25
                    
                    line = fo.readline()
                
                motifs.append(Motif(pfm[:]))
                motifs[-1].id = id
                motifs[-1].align = align[:]
            line = fo.readline()

        return motifs

__tools__ = {
        "xxmotif": XXmotif,
        "homer": Homer, 
        "bioprospector":BioProspector,
        "hms": Hms,
        "amd": Amd,
        "improbizer": Improbizer,
        "trawler": Trawler,
        "weeder": Weeder,
        "motifsampler": MotifSampler,
        "mdmodule": MDmodule,
        "chipmunk": ChIPMunk,
        "posmo": Posmo,
        "gadem": Gadem,
        "jaspar": Jaspar,
        "meme": Meme,
        "memew": MemeW,
    }

