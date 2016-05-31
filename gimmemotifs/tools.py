# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Interface module for all motif programs """

# Python imports
import re
import os
import sys
from subprocess import Popen, PIPE, call
import shutil
from tempfile import NamedTemporaryFile, mkdtemp
import StringIO

# gimme imports
from gimmemotifs.config import MotifConfig
from gimmemotifs.fasta import Fasta
from gimmemotifs.utils import which

# Necessary for when setup.py needs to import 
# refactor at some point
try:
    from gimmemotifs.motif import read_motifs, Motif
except ImportError:
    pass


def locate_tool(tool, verbose=True): 
    tool = re.sub(r'[^a-zA-Z]','',tool) 
    m = eval(tool)() 
    tool_bin = which(m.cmd) 
    if tool_bin:
        if verbose:
            print "Found {} in {}".format(m.name, tool_bin) 
        return tool_bin 
    else: 
        print "Couldn't find {}".format(m.name)

class MotifProgram(object):
    config = MotifConfig()
    local_bin = None

    def __init__(self):
        pass

    def bin(self):
        if self.local_bin:
            return self.local_bin
        else:
            return self.config.bin(self.name)

    def dir(self):
        return self.config.dir(self.name)

    def is_configured(self):
        return self.config.is_configured(self.name)
    
    def is_installed(self):
        return self.is_configured() and os.access(self.bin(), os.X_OK)

    def run(self, fastafile, savedir, params=None, tmp=None):

        if not self.is_configured():
            raise ValueError("%s is not configured" % self.name)

        if not self.is_installed():
            raise ValueError("%s is not installed or not correctly configured" % self.name)
        
        self.tmpdir = mkdtemp(prefix="{0}.".format(self.name), dir=tmp)
 
        try:
            return self._run_program(self.bin(), fastafile, savedir, params)
        except KeyboardInterrupt:
            return ([], "Killed", "Killed")
#        except Exception as e:
#            return ([], "", e.strerror)

class XXmotif(MotifProgram):

    def __init__(self):
        self.name = "XXmotif"
        self.cmd = "XXmotif"
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir="", params=None):
        if params is None:
            params = {}
        
        default_params = {"single":False, "background":None, "analysis":"medium", "number":5, "width":10}
        default_params.update(params)
        
        cmd = bin
        
        fastafile = os.path.abspath(fastafile)
        
        bgfile = os.path.abspath(default_params["background"])
        background = ""
        if bgfile:
            background = " --negSet {0} ".format(bgfile)

        outfile = os.path.join(self.tmpdir, os.path.basename(fastafile.replace(".fa", ".pwm")))
        
        stdout = ""
        stderr = ""
        
        strand = ""
        if not default_params["single"]:
            strand = " --revcomp "

        cmd = "%s %s %s --localization --batch --no-graphics %s %s" % (
            cmd,
            self.tmpdir, 
            fastafile,
            background,
            strand
            )

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        motifs = []
        
        if os.path.exists(outfile):
            motifs = read_motifs(open(outfile), fmt="xxmotifs")
            for m in motifs:
                m.id = "{0}_{1}".format(self.name, m.id)
        else:
            stdout += "\nMotif file {0} not found!\n".format(outfile)
            stderr += "\nMotif file {0} not found!\n".format(outfile)
        
        return motifs, stdout, stderr
        
class Homer(MotifProgram):

    def __init__(self):
        self.name = "Homer"
        self.cmd = "homer2"
        self.use_width = True

    def _run_program(self, bin, fastafile, savedir="", params=None):
        
        default_params = {"single":False, "background":None, "analysis":"medium", "number":5, "width":10}
        if params is not None: 
            default_params.update(params)
        
        homer = bin
        
        fastafile = os.path.abspath(fastafile)
        
        # Background file is essential!
        if not default_params["background"]:
            print "Background file needed!"
            sys.exit()
        
        bgfile = os.path.abspath(default_params["background"])
        
        outfile = NamedTemporaryFile(
                dir=self.tmpdir, 
                prefix= "homer_w{}.".format(default_params["width"])
                ).name
        
        stderr = ""
        
        strand = ""
        if default_params["single"]:
            strand = " -strand + "

        cmd = "%s denovo -i %s -b %s -len %s -S %s %s -o %s -p 8" % (
            homer,
            fastafile,
            bgfile,
            default_params["width"],
            default_params["number"],
            strand,
            outfile)

        stdout = "Running command:\n{}\n".format(cmd)
        
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=self.tmpdir) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        motifs = []
        
        if os.path.exists(outfile):
            motifs = read_motifs(open(outfile), fmt="pwm")
            for i, m in enumerate(motifs):
                m.id = "{}_{}_{}".format(self.name, default_params["width"], i + 1)
        
        return motifs, stdout, stderr

class BioProspector(MotifProgram):
    def __init__(self):
        self.name = "BioProspector"
        self.cmd = "BioProspector"
        self.use_width = True

    def _run_program(self, bin, fastafile, savedir="", params=None):
        
        default_params = {"single":False, "background":None, "analysis":"medium", "number":5, "width":10}
        if params is not None: 
            default_params.update(params)
        
        prospector = bin
        
        fastafile = os.path.abspath(fastafile)
        
        # Background file is essential!
        if not default_params["background"]:
            print "Background file needed!"
            sys.exit()
        
        #bgfile = os.path.abspath(default_params["background"])
        outfile = os.path.join(self.tmpdir, "bioprospector.out")    
        
        stdout = ""
        stderr = ""
        
        strand = 2
        if default_params["single"]:
            strand = 1

        cmd = "%s -i %s -W %s -d %s -b %s -r %s -o %s" % (
            prospector,
            fastafile,
            default_params["width"],
            strand,
            default_params["background"],
            default_params["number"],
            outfile)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        motifs = []
        
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()
        
        return motifs, stdout, stderr

    def parse(self, fo):
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
    def __init__(self):
        self.name = "HMS"
        self.cmd = "hms"
        self.use_width = False
    
    def _run_program(self, bin, fastafile, savedir="", params=None):

        hms = bin
        thetas = ["theta%s.txt" % i for i in [0,1,2,3]]

        fastafile = os.path.abspath(fastafile)
        
        fgfile = os.path.join(self.tmpdir, "HMS.in.fa")
        summitfile = os.path.join(self.tmpdir, "HMS.in.summits.txt")
        outfile = os.path.join(self.tmpdir, "thetafinal.txt")    
    
        hmsdir = os.path.join(self.config.get_tools_dir(), "HMS")
        shutil.copy(fastafile, fgfile)
        for t in thetas:
            shutil.copy(os.path.join(hmsdir, t), self.tmpdir)
        
        fa = Fasta(fgfile)
        out = open(summitfile, "w")
        for seq in fa.seqs:
            out.write("%s\n" % (len(seq) / 2))
        out.close()
        
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        stdout = ""
        stderr = ""
    
        cmd = "%s -i %s -w 21 -dna 4 -iteration 50 -chain 20 -seqprop -0.1 -strand 2 -peaklocation %s -t_dof 3 -dep 2" % (hms, fgfile, summitfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        os.chdir(current_path)
        motifs = []
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()
        
        return motifs, stdout, stderr

    def parse(self, fo):
        motifs = []
        m = [[float(x) for x in fo.readline().strip().split(" ")] for i in range(4)]
        matrix = [[m[0][i], m[1][i],m[2][i],m[3][i]] for i in range(len(m[0]))]
        motifs = [Motif(matrix)]
        motifs[-1].id = self.name
        
        return motifs

class Amd(MotifProgram):
    def __init__(self):
        self.name = "AMD"
        self.cmd = "AMD.bin"
        self.use_width = False
    
    def _run_program(self, bin, fastafile, savedir="", params=None):
        
        default_params = {"background":None}
        if params is not None: 
            default_params.update(params)
        
        amd = bin
        
        fastafile = os.path.abspath(fastafile)
        
        # Background file is essential!
        if not default_params["background"]:
            print "Background file needed!"
            sys.exit(1)
        
        bgfile = os.path.abspath(default_params["background"])
        fgfile = os.path.join(self.tmpdir, "AMD.in.fa")
        outfile = fgfile + ".Matrix"    
    
        shutil.copy(fastafile, fgfile)
        print fgfile
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        stdout = ""
        stderr = ""
    
        cmd = "%s -F %s -B %s" % (amd, fgfile, bgfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        os.chdir(current_path)
        motifs = []
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()
        
        return motifs, stdout, stderr

    def parse(self, fo):
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
    def __init__(self):
        self.name = "Improbizer"
        self.cmd = "ameme"
        self.use_width = False 

    def _run_program(self, bin, fastafile, savedir="", params=None):
        
        default_params = {"background":None, "number":10}
        if params is not None: 
            default_params.update(params)
        
        ameme = bin
        
        fastafile = os.path.abspath(fastafile)
        
        # Not strictly necessary, but recommended
        if not default_params["background"]:
            print "Background file needed!"
            sys.exit()
        
        bgfile = os.path.abspath(default_params["background"])
        outfile = os.path.join(self.tmpdir, "improbizer.out.html")    
        
        current_path = os.getcwd()
        os.chdir(self.tmpdir)
        
        stdout = ""
        stderr = ""
        cmd = "%s good=%s bad=%s numMotifs=%s > %s" % (ameme, fastafile, bgfile, default_params["number"], outfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        motifs = []
        if os.path.exists(outfile):
            f = open(outfile)
            motifs = self.parse(f)
            f.close()
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
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
    def __init__(self):
        self.name = "trawler"
        self.cmd = "trawler.pl"
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir="", params=None):
        
        default_params = {"single":False, "background":None}
        if params is not None: 
            default_params.update(params)
        
        trawler = bin
        
        fastafile = os.path.abspath(fastafile)
        if not default_params["background"]:
            print "Background file needed!"
            sys.exit()
        bgfile = os.path.abspath(default_params["background"])
        savedir = os.path.abspath(savedir)
        
        #savedir = "/tmp/trawler/"

        tmp = NamedTemporaryFile(dir=self.tmpdir, delete=False)
        shutil.copy(fastafile, tmp.name)
        fastafile = tmp.name
    
        current_path = os.getcwd()
        os.chdir(self.dir())
        
        stdout = ""
        stderr = ""
        strand = "double"
        if default_params["single"]:
            strand = "single"
        cmd = "%s -sample %s -background %s -directory %s -strand %s" % (trawler, fastafile, bgfile, self.tmpdir, strand)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        os.chdir(current_path)
        motifs = []
        out_name = [dir for dir in os.listdir(self.tmpdir) if dir.startswith("tmp")][-1]
        out_file = os.path.join(self.tmpdir, out_name, "result", "%s.pwm" % out_name)
        if os.path.exists(out_file):
            motifs = read_motifs(open(os.path.join(
                                                    self.tmpdir, 
                                                    out_name, 
                                                    "result", 
                                                    "%s.pwm" % out_name)),
                                                    fmt="pwm")
        
        # remove temporary files
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)
        
        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        return []

def run_weeder_subset(weeder, fastafile, w, e, organism, strand):
    cmd = "%s -f %s -W %s -e %s -R 50 -O %s %s" % (weeder, fastafile, w, e, organism, strand)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
    out,err = p.communicate()
    return out, err

class Weeder(MotifProgram):
    def __init__(self):
        self.name = "Weeder"
        self.cmd = "weederTFBS.out"
        self.use_width = False

    def _run_program(self, bin,fastafile, savedir="", params=None):
        #try: 
        #    from gimmemotifs.mp import pool
        #except:
        #    pass


        default_params = {"analysis":"small", "organism":"hg18", "single":False, "parallel":True}
        if params is not None: 
            default_params.update(params)
        
        organism = default_params["organism"]
        weeder_organism = ""
        weeder_organisms = {
            "hg18":"HS",
            "hg19":"HS",
            "mm9":"MM",
            "rn4":"RN",
            "dm3":"DM",
            "fr2": "FR",
            "danRer6": "DR",
            "danRer7": "DR",
            "galGal3": "GG",
            "ce3": "CE",
            "anoGam1": "AG",
            "yeast":"SC",
            "sacCer2":"SC",
            "xenTro2":"XT",
            "xenTro3":"XT"}
        if weeder_organisms.has_key(organism):
            weeder_organism = weeder_organisms[organism]
        else:
            return []    

        weeder = bin
        adviser = weeder.replace("weederTFBS", "adviser")
    
        
        weeder_dir = bin.replace("weederTFBS.out", "")
        if self.is_configured():
            weeder_dir = self.dir()

        freq_files = os.path.join(weeder_dir, "FreqFiles")
        if not os.path.exists(freq_files):
            raise ValueError, "Can't find FreqFiles directory for Weeder"
                

        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)

        tmp = NamedTemporaryFile(dir=self.tmpdir)
        name = tmp.name
        tmp.close()
        shutil.copy(fastafile, name)
        fastafile = name
    
        current_path = os.getcwd()
        os.chdir(weeder_dir)
        
        coms = ((8,2),(6,1))

        strand = "-S"
        if default_params["single"]:
            strand = ""
            
        if default_params["analysis"] == "xl":
             coms = ((12,4),(10,3),(8,2),(6,1))
        elif default_params["analysis"] == "large":
             coms = ((10,3),(8,2),(6,1))
        elif default_params["analysis"] == "medium":
             coms = ((10,3),(8,2),(6,1))
        
        # TODO: test organism
        stdout = ""
        stderr = ""
        
        default_params["parallel"] = False 
        if default_params["parallel"]:
            jobs = []
            #for (w,e) in coms:
            #    jobs.append(pool.apply_async(
            #        run_weeder_subset, 
            #        (weeder, fastafile, w, e, weeder_organism, strand,)
            #        ))

            #for job in jobs:
            #    out,err = job.get()
            #    stdout += out
            #    stderr += err
        else:

            for (w,e) in coms:
                out,err = run_weeder_subset(weeder, fastafile, w, e, weeder_organism, strand)
                stdout += out
                stderr += err
    
        cmd = "%s %s" % (adviser, fastafile)
        #print cmd
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        out,err = p.communicate()
        stdout += out
        stderr += err
        
        os.chdir(current_path)

        motifs = []
        if os.path.exists(fastafile + ".wee"):
            f = open(fastafile + ".wee")
            motifs = self.parse(f)
            f.close()
        
        for ext in [".wee", ".html", ".mix", ""]:
            if os.path.exists(fastafile + ext):
                os.unlink(fastafile + ext)

        return motifs, stdout, stderr

    def parse(self, fo):
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

        p = re.compile(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
        pa = re.compile(r'(\s*\d+\s+.\s+([ACGT]+)\s+.+\))')
        pwm = []
        align = []
        c = 1
        for line in fo.readlines():
            m = p.search(line)
            if m:
                pwm.append([int(m.group(x)) for x in [2,3,4,5]])
            
            m = pa.search(line)
            if m:
                align.append(m.group(2))

            elif line.startswith("===="):
                motifs.append(Motif())
                #total = float(pwm[0][0] + pwm[0][1] + pwm[0][2] + pwm[0][3])
                #motifs[-1].pwm = [[x / total for x in row] for row in pwm]
                motifs[-1].id = "Weeder_%s" % c
                motifs[-1].align= align[:]
                
                width = len(align[0])
                pfm =  [[0 for x in range(4)] for x in range(width)]
                for row in align:
                    for i in range(len(row)):
                        pfm[i][nucs[row[i]]] += 1
                total = float(len(align)) 
                pwm = [[((x)/total) for x in row] for row in pfm]
                motifs[-1].pwm = pwm[:]
                motifs[-1].pfm = pfm[:]
                
                align = []
                c += 1
                pwm = []
                    
        return motifs

class MotifSampler(MotifProgram):
    def __init__(self):
        self.name = "MotifSampler"
        self.cmd = "MotifSampler"
        self.use_width = True

    def _run_program(self, bin, fastafile, savedir, params=None):
        
        default_params = {"width":10, "background":"", "single":False, "number":10}
        if params is not None: 
            default_params.update(params)
        
        background = default_params['background']
        width = default_params['width']
        number = default_params['number']

        if not background:
            if default_params["organism"]:
                org = default_params["organism"]
                background = os.path.join(
                        self.config.get_bg_dir(), 
                        "{}.{}.bg".format(org, "MotifSampler"))
            else:            
                raise Exception, "No background specified for {}".format(self.name)

        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
    
        tmp = NamedTemporaryFile(dir=self.tmpdir)
        pwmfile = tmp.name

        tmp2  = NamedTemporaryFile(dir=self.tmpdir)
        outfile = tmp2.name
    
        strand = 1
        if default_params["single"]:
            strand = 0

        # TODO: test organism
        cmd = "%s -f %s -b %s -m %s -w %s -n %s -o %s -s %s > /dev/null 2>&1" % (bin, fastafile, background, pwmfile, width, number, outfile, strand)
        #print cmd
        #p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        #stdout, stderr = p.communicate()
        
        stdout,stderr = "",""
        p = Popen(cmd, shell=True)
        p.wait()

        motifs = []
        #if os.path.exists(pwmfile):
        #    motifs = self.parse(open(pwmfile))
        if os.path.exists(outfile):
            motifs = self.parse_out(open(outfile))
        
        # remove temporary files
        tmp.close()
        tmp2.close()
        
        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)
        
        return motifs, stdout, stderr

    def parse(self, fo):
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
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}
        pseudo = 0.0 # Should be 1/sqrt(# of seqs)
        aligns = {}
        for line in fo.readlines():
            if line.startswith("#"):
                pass
            elif len(line) > 1:
                vals = line.strip().split("\t")
                id, site = [x.strip().split(" ")[1].replace('"',"") for x in vals[8].split(";") if x]
                #if vals[6] == "+":
                if site.upper().find("N") == -1:
                    aligns.setdefault(id, []).append(site)
                #else:
                #    print site, rc(site)
                #    aligns.setdefault(id, []).append(rc(site))
                        
        for id, align in aligns.items():
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
            m.id = id
            motifs.append(m)    
        return motifs

class MDmodule(MotifProgram):
    def __init__(self):
        self.name = "MDmodule"
        self.cmd = "MDmodule"
        self.use_width = True
        
    def _run_program(self, bin, fastafile, savedir, params=None):
        
        default_params = {"width":10, "number":10}
        if params is not None: 
            default_params.update(params)
        
        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
    
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
        
        stdout = "cmd: {}\n".format(cmd) + stdout 
            
        motifs = []
        if os.path.exists(pwmfile):
            motifs = self.parse(open(pwmfile))
        
        os.chdir(current_path)
        
        for motif in motifs:
            motif.id = "%s_%s" % (self.name, motif.id)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}
        p = re.compile(r'(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
        pf = re.compile(r'>.+\s+[bf]\d+\s+(\w+)')

        pwm = []
        pfm = []
        align = []
        id = ""
        for line in fo.readlines():
            if line.startswith("Motif"):
                if id:
                    motifs.append(Motif())
                    motifs[-1].id = id
                    motifs[-1].pwm = pwm
                    motifs[-1].pfm = pfm
                    motifs[-1].align = align
                    pwm = []
                    pfm = []
                    align = []
                id = line.split("\t")[0]
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
            motifs[-1].id = id
            motifs[-1].pwm = pwm
            motifs[-1].pfm = pfm
            motifs[-1].align = align

        return motifs

class ChIPMunk(MotifProgram):
    def __init__(self):
        self.name = "ChIPMunk"
        self.cmd = "ChIPMunk.sh"
        self.use_width = True

    def _run_program(self, bin, fastafile, savedir, params=None):

        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)

        basename = "munk_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        out = open(new_file, "w")
        f = Fasta(fastafile)
        for name,seq in f.items():
            header = " ".join(["%0.1f" % x for x in range(len(seq) / 2) + range(len(seq) / 2, 0, -1)])
            out.write(">%s\n" % header)
            out.write("%s\n" % seq)
        out.close()
        
        fastafile = new_file
        outfile = fastafile + ".out"

        current_path = os.getcwd()
        os.chdir(self.dir())
        
        cmd = "%s %s %s yes 1.0 p:%s > %s" % (bin, params["width"], params["width"], fastafile, outfile)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
        stdout, stderr = p.communicate()

        motifs = []
        if os.path.exists(outfile):
            motifs = self.parse(open(outfile))
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr
        
    def parse(self, fo):
        #KDIC|6.124756232026243
        #A|517.9999999999999 42.99999999999999 345.99999999999994 25.999999999999996 602.9999999999999 155.99999999999997 2.9999999999999996 91.99999999999999
        #C|5.999999999999999 4.999999999999999 2.9999999999999996 956.9999999999999 91.99999999999999 17.999999999999996 22.999999999999996 275.99999999999994
        #G|340.99999999999994 943.9999999999999 630.9999999999999 6.999999999999999 16.999999999999996 48.99999999999999 960.9999999999999 14.999999999999998
        #T|134.99999999999997 7.999999999999999 19.999999999999996 9.999999999999998 287.99999999999994 776.9999999999999 12.999999999999998 616.9999999999999
        #N|999.9999999999998
        line = fo.readline()
        while not line.startswith("A|"):
            line = fo.readline() 
        matrix = []
        for i in range(4):
            matrix.append([float(x) for x in line.strip().split("|")[1].split(" ")])
            line = fo.readline()
        #print matrix
        matrix = [[matrix[x][y] for x in range(4)] for y in range(len(matrix[0]))]
        #print matrix
        m = Motif(matrix)
        m.id = "ChIPMunk_w%s" % len(m)
        return [m]


class Posmo(MotifProgram):
    def __init__(self):
        self.name = "Posmo"
        self.cmd = "posmo"
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir, params=None):
        
        default_params = {}
        if params is not None: 
            default_params.update(params)
        
        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
    
        basename = "posmo_in.fa"

        new_file = os.path.join(self.tmpdir, basename)
        shutil.copy(fastafile, new_file)
        
        fastafile = new_file
        pwmfile = fastafile + ".pwm"
    
        motifs = []
        current_path = os.getcwd()
        os.chdir(self.tmpdir)    
        for x in ["111111", "11111111"]:
            outfile = "%s.%s.out" % (fastafile, x)
            cmd = "%s 5000 %s %s 1.6 2.5 20 200" % (bin, x, fastafile)
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
            stdout, stderr = p.communicate()
    
            context_file = fastafile.replace(basename, "context.%s.%s.txt" % (basename, x))
            cmd = "%s %s %s simi.txt 0.88 10 2 10" % (bin.replace("posmo","clusterwd"), context_file, outfile)
            p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
            out, err = p.communicate()
            stdout += out
            stderr += err
        
            if os.path.exists(outfile):
                motifs += self.parse(open(outfile))
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

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
    def __init__(self):
        self.name = "GADEM"
        self.cmd = "gadem"
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir, params=None):
        
        default_params = {}
        if params is not None: 
            default_params.update(params)
        
        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
    
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
            motifs = self.parse(open(pwmfile))
        
        os.chdir(current_path)
        
        return motifs, stdout, stderr

    def parse(self, fo):
        motifs = []
        nucs = {"A":0,"C":1,"G":2,"T":3}

        lines = fo.readlines()
        for i in range(0, len(lines), 5):
            align = []
            pwm = []
            pfm = []
            id = ""
            line = lines[i].strip()
            id = line[1:]
            number = id.split("_")[0][1:]
            if os.path.exists("%s.seq" % number):
                for l in open("%s.seq" % number).readlines():
                    if not "x" in l and not "n" in l:
                        l = l.strip().upper()
                        align.append(l)
                        if not pfm:
                            pfm = [[0 for x in range(4)] for x in range(len(l))]
                        for p in range(len(l)):
                            pfm[p][nucs[l[p]]] += 1

            m = [l.strip().split(" ")[1].split("\t") for l in lines[i + 1: i + 5]]

            pwm = [[float(m[x][y]) for x in range(4)] for y in range(len(m[0]))]


            motifs.append(Motif(pwm))
            motifs[-1].id = id
            #motifs[-1].pwm = pwm
            if align:
                pass
                motifs[-1].pfm = pfm
                motifs[-1].align = align

        return motifs

class Jaspar(MotifProgram):
    def __init__(self):
        self.name = "JASPAR"
        self.cmd = "/bin/false"    
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir, params=None):
        fname = os.path.join(self.config.get_motif_dir(), "JASPAR2010_vertebrate.pwm")
        motifs =  read_motifs(open(fname), fmt="pwm")
        for motif in motifs:
            motif.id = "JASPAR_%s" % motif.id
        return motifs, "", ""

class Meme(MotifProgram):
    def __init__(self):
        self.name = "MEME"
        self.cmd = "meme.bin"
        self.use_width = True

    def _run_program(self, bin, fastafile, savedir, params=None):
        #EVT = 1.0
        default_params = {"width":10, "single":False, "number":10}
        if params is not None: 
            default_params.update(params)
        
        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
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
        motifs = self.parse(StringIO.StringIO(stdout))
        
        # Delete temporary files
        tmp.close()
         
        return motifs, stdout, stderr

    def parse(self, fo):
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
                            if l[pos] in nucs.keys():
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
    def __init__(self):
        self.name = "MEMEW"
        self.cmd = "meme.bin"
        self.use_width = False

    def _run_program(self, bin, fastafile, savedir, params=None):
        
        #EVT = 1.0
        default_params = {"single":False, "number":10}
        if params is not None: 
            default_params.update(params)
        
        fastafile = os.path.abspath(fastafile)
        savedir = os.path.abspath(savedir)
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
        motifs = self.parse(StringIO.StringIO(stdout))
        
        # Delete temporary files
        tmp.close()
         
        return motifs, stdout, stderr

    def parse(self, fo):
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
                            if l[pos] in nucs.keys():
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
