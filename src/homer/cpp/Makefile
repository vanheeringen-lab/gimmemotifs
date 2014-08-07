DEBUG = -O3 -Wall
#DEBUG = -O3 -Wall
#DEBUG = -O3 -march=pentium4 -fomit-frame-pointer
#DEBUG = -O3 -D_FILE_OFFSET_BITS=64

#for linux
COMPILER = g++
LIBRARIES = -lpthread -lstdc++ -lm
#for SunOS
#COMPILER = CC
#LIBRARIES = -lCstd -lm

homer2 homer getMappableRegions makeTagDirectory getPeakTags makeUCSCfile getDifferentialPeaks mergePeaks genomeOntology homerTools analyzeHiC getGenomeTilingPeaks : Hashtable.o Motif.o statistics.o CommandLine.o homer.o getMappableRegions.o SeqTag.o makeTagDirectory.o getPeakTags.o findPeaks.o makeUCSCfile.o getDifferentialPeaks.o assignGenomeAnnotation.o mergePeaks.o genomeOntology.o homerTools.o analyzeHiC.o Motif2.o homer2.o Clustering.o getGenomeCorrelation.o
	$(COMPILER) homer.o CommandLine.o Motif.o statistics.o Hashtable.o $(LIBRARIES) $(DEBUG) -o homer
	$(COMPILER) getMappableRegions.o Hashtable.o $(LIBRARIES) $(DEBUG) -o getMappableRegions
	$(COMPILER) getPeakTags.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o getPeakTags
	$(COMPILER) makeTagDirectory.o SeqTag.o Motif2.o Clustering.o  Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o makeTagDirectory
	$(COMPILER) findPeaks.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o findPeaks
	$(COMPILER) makeUCSCfile.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o makeUCSCfile
	$(COMPILER) getDifferentialPeaks.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o getDifferentialPeaks
	$(COMPILER) assignGenomeAnnotation.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o assignGenomeAnnotation
	$(COMPILER) mergePeaks.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o mergePeaks
	$(COMPILER) genomeOntology.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o genomeOntology
	$(COMPILER) homerTools.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o homerTools
	$(COMPILER) analyzeHiC.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o analyzeHiC
	$(COMPILER) getGenomeCorrelation.o SeqTag.o Motif2.o Clustering.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o getGenomeTilingPeaks
	$(COMPILER) homer2.o Motif2.o Hashtable.o statistics.o $(LIBRARIES) $(DEBUG) -o homer2
	cp findPeaks ../bin/
	cp mergePeaks ../bin/
	cp genomeOntology ../bin/
	cp getPeakTags ../bin/
	cp makeUCSCfile ../bin/
	cp getDifferentialPeaks ../bin/
	cp assignGenomeAnnotation ../bin/
	cp homerTools ../bin/
	cp analyzeHiC ../bin/
	cp homer ../bin/
	cp makeTagDirectory ../bin/
	cp homer2 ../bin/
	cp getMappableRegions ../bin/
	cp getGenomeTilingPeaks ../bin/getGenomeTilingPeaks

homerTools.o : homerTools.cpp homerTools.h
	$(COMPILER) homerTools.cpp -c $(DEBUG)
getGenomeCorrelation.o : getGenomeCorrelation.cpp SeqTag.h
	$(COMPILER) getGenomeCorrelation.cpp -c $(DEBUG)
Clustering.o : Clustering.cpp Clustering.h
	$(COMPILER) Clustering.cpp -c $(DEBUG)
analyzeHiC.o : analyzeHiC.cpp
	$(COMPILER) analyzeHiC.cpp -c $(DEBUG)
Motif2.o : Motif2.cpp Motif2.h
	$(COMPILER) Motif2.cpp -c $(DEBUG)
getMappableRegions.o : getMappableRegions.cpp
	$(COMPILER) getMappableRegions.cpp -c $(DEBUG)
Hashtable.o : Hashtable.cpp Hashtable.h
	$(COMPILER) Hashtable.cpp -c $(DEBUG)
statistics.o : statistics.cpp statistics.h
	$(COMPILER) statistics.cpp -c $(DEBUG) 
Motif.o : Motif.cpp Motif.h
	$(COMPILER) Motif.cpp -c $(DEBUG)
CommandLine.o : CommandLine.cpp CommandLine.h
	$(COMPILER) CommandLine.cpp -c $(DEBUG)
homer.o : homer.cpp
	$(COMPILER) homer.cpp -c $(DEBUG)
homer2.o : homer2.cpp
	$(COMPILER) homer2.cpp -c $(DEBUG)
getPeakTags.o : getPeakTags.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) getPeakTags.cpp -c $(DEBUG)
makeUCSCfile.o : makeUCSCfile.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) makeUCSCfile.cpp -c $(DEBUG)
getDifferentialPeaks.o : getDifferentialPeaks.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) getDifferentialPeaks.cpp -c $(DEBUG)
mergePeaks.o : mergePeaks.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) mergePeaks.cpp -c $(DEBUG)
genomeOntology.o : genomeOntology.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) genomeOntology.cpp -c $(DEBUG)
findPeaks.o : findPeaks.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) findPeaks.cpp -c $(DEBUG)
makeTagDirectory.o : makeTagDirectory.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) makeTagDirectory.cpp -c $(DEBUG)
assignGenomeAnnotation.o : assignGenomeAnnotation.cpp SeqTag.h SeqTag.cpp
	$(COMPILER) assignGenomeAnnotation.cpp -c $(DEBUG)
SeqTag.o : SeqTag.cpp SeqTag.h
	$(COMPILER) SeqTag.cpp -c $(DEBUG)

clean:
	rm -f *.o;
	rm -f getMappableRegions makeTagDirectory findPeaks mergePeaks genomeOntology getPeakTags makeUCSCfile getDifferentialPeaks assignGenomeAnnotation homerTools homer homer2 analyzeHiC getGenomeTilingPeaks;
