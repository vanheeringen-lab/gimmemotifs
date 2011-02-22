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
$location = "http://ws.nbcr.net/opal2/services/MEME_4.4.0";
$fastaFile = "./crp0.fasta";

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
# Use DNA alphabet, and find 3 motifs.
$commandline = "meme crp0.fasta -dna -nmotifs 3";
$req->setArgs($commandline);

# Set name and content of input file.
# Only one input file is used, the FASTA file.
open MYFILE, $fastaFile;
my $file = do {local $/; <MYFILE>};
$inputFile = InputFileType->new($fastaFile, $file);

# Attache input file to job.
$req->setInputFile($inputFile);

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
my $output_dir = "meme_out";
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
  my $name = $_->getName();
  if ($name eq "bin" or $name eq "data") {
    # Don't try to retreive "bin", or "data", there directories.
    next;
  }
  print "\t", $name, ": ", $_->getURL(),"\n";
  $req = OutputsByNameInputType->new($jobid, $name);
  $result = $memeService->getOutputAsBase64ByName($req);
  $outfile = $output_dir . "/" . $name;
  open OUTFILE, "> $outfile";
  print OUTFILE $result;
  close OUTFILE;
}

exit 0;
