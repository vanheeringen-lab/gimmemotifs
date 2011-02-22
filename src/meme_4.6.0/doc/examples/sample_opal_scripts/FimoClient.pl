# This script submits a FASTA file to the MEME web service at NBCR
# and downloads the resulting output to the directory meme_out.
# Based on the script pdb2pqrclient.pl included in the 1.92 Opal-Perl
# toolkit, available at http://www.nbcr.net/software/opal/

# Import the Opal libraries
use OpalServices;
use OpalTypes;
use File::Basename;

# Set the location of our service as well as the name of the FASTA file to input
# A list of NBCR services can be found at http://ws.nbcr.net/opal2/dashboard?command=serviceList
$location = "http://ws.nbcr.net/opal2/services/FIMO_4.4.0";

# Instantiate a new service object to interact with.  Pass it the service location
$memeService = OpalServices->new(service_url => $location);

# Make a remote invocation using the service object to retrieve the application metadata
# including the usage statement.
$result = $memeService->getAppMetadata();
print "App Metadata:\n";
print "\tUsage: ", $result->getUsage(),"\n";
print "\tInfo:\n";
print $result->getInfo(), "\n";

# Similar invocation for getting the application configuration
$result = $memeService->getAppConfig();
print "\nApp Config:\n";
print "\tBinary Location: ", $result->getBinaryLocation(),"\n";
print "\tDefault Args: ", $result->getDefaultArgs(),"\n";

# Instantiate a new job request
$req = JobInputType->new();

# Set the command-line arguments for the job.
# The first two items must be the motif filename
# and the sequence db filename. Note that the 
# sequence db filename must be set to uploaded_db
# The FIMO web service uses this to indicate that
# The sequence db is being uploaded, rather then
# using one of the preloaded databases.
# the '--oc .' option is required for Opal
# to see all the output files created by FIMO
$commandline = "crp0.meme.xml crp0.fasta --oc .";
# Alternatively, use the argList setting below to scan the 
# preloaded Saccaromyces cerevisiae genome.
#$commandline = "crp0.meme.xml saccharomyces_cerevisiae.na  --oc ."
$req->setArgs($commandline);

# Get contents of local files to be used as input
$fastaFile = "./crp0.fasta";
$motifFile = "./crp0.meme.xml";
open MYFILE, $fastaFile;
my $fasta = do {local $/; <MYFILE>};
close MYFILE;
open MYFILE, $motifFile;
my $motif = do {local $/; <MYFILE>};
close MYFILE;

# Set name and content of remote input file.
# Two input files are used, the motif file
#and the FASTA file.
$fastaInputFile = InputFileType->new($fastaFile, $fasta);
$motifInputFile = InputFileType->new($motifFile, $motif);
$req->setInputFile($motifInputFile, $fastaInputFile);

# Launch the job and retrieve job ID
$result = $memeService->launchJob($req);
print "Job Launched:\n";
print "\tCode: ", $result->getCode(),"\n";
print "\tMessage: ", $result->getMessage(),"\n";
print "\tBase URL: ", $result->getBaseURL(),"\n";

# Loop until job is complete (GRAM status 8 = Finished)
$statuscode = 0;
$jobid = $result->getJobID();
while ($statuscode != 8) {
  $status = $memeService->queryStatus($jobid);
  $statuscode=$status->getCode();
  print "Query Status:\n";
  print "\tCode: ", $status->getCode(),"\n";
  print "\tMessage: ", $status->getMessage(),"\n";
  print "\tBase URL: ", $status->getBaseURL(),"\n";
  sleep(5);
}

# Download the output files for the job.

# Create output directory
my $output_dir = "fimo_out";
mkdir $output_dir;

# Retreive output object
$output = $memeService->getOutputs($jobid);

# The MEME services writes stdout and stderr to
# files that don't show up in the list of output
# files, so we handle stdout and stderr separately.
print "Downloading standard out\n";
$req = OutputsByNameInputType->new($jobid, "stdout.txt");
$result = $memeService->getOutputAsBase64ByName($req);
$outfile = $output_dir . "/stdout.txt";
open OUTFILE, "> $outfile";
print OUTFILE $result;
close OUTFILE;
print "Downloading standard error\n";
$req = OutputsByNameInputType->new($jobid, "stderr.txt");
$result = $memeService->getOutputAsBase64ByName($req);
$outfile = $output_dir . "/stderr.txt";
open OUTFILE, "> $outfile";
print OUTFILE $result;
close OUTFILE;

# Download all other output files
print "\nDownloading output files - : \n";
@list = $output->getFiles();
foreach (@list) {
  if ($_->getName() eq "bin") {
    # Don't try to retreive "bin", it's a directory.
    next;
  }
  print "\t", $_->getName(), ": ", $_->getURL(),"\n";
  $req = OutputsByNameInputType->new($jobid, $_->getName());
  $result = $memeService->getOutputAsBase64ByName($req);
  $outfile = $output_dir . "/" . $_->getName();
  open OUTFILE, "> $outfile";
  print OUTFILE $result;
  close OUTFILE;
}

exit 0;
