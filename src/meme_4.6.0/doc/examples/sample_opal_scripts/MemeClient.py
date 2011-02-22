# This script submits a FASTA file to the MEME web service at NBCR
# and downloads the resulting output to the directory meme_out.
# Based on the script pdb2pqrclient.pl included in the 2.0.0 Opal-Python
# toolkit, available at http://www.nbcr.net/software/opal/

# Import the Opal libraries
from AppService_client import \
     AppServiceLocator, getAppMetadataRequest, launchJobRequest, \
     queryStatusRequest, getOutputsRequest, \
     launchJobBlockingRequest, getOutputAsBase64ByNameRequest
from AppService_types import ns0
from os import mkdir
from time import sleep
from ZSI.TC import String

# Set the location of our service
# A list of NBCR services can be found at http://ws.nbcr.net/opal2/dashboard?command=serviceList
serviceURL = "http://ws.nbcr.net/opal2/services/MEME_4.4.0"

# Instantiate a new service object to interact with.  Pass it the service location
appLocator = AppServiceLocator()
appServicePort = appLocator.getAppServicePort(serviceURL)
    
# Instantiate a new non-blocking job request
req = launchJobRequest()

# Set the command-line arguments for the job.
# Use DNA alphabet, and find 3 motifs.
req._argList = "meme crp0.fasta -dna -mod zoops -nmotifs 3"

# Get contents of local file to be used as input
fastaFile = open("./crp0.fasta", "r")
fastaFileContents = fastaFile.read()
fastaFile.close()

# Set name and content of remote input file.
# Only one input file is used, the FASTA file.
inputFiles = []
inputFile = ns0.InputFileType_Def('inputFile')
inputFile._name = 'crp0.fasta'
inputFile._contents = fastaFileContents
inputFiles.append(inputFile)
req._inputFile = inputFiles

# Launch non-blocking job and retrieve job ID
print "Launching remote MEME job"
resp = appServicePort.launchJob(req)
jobID = resp._jobID
print "Received Job ID:", jobID

# Poll for job status
status = resp._status
print "Polling job status"
while 1:
  # print current status
  print "Status:"
  print "\tCode:", status._code
  print "\tMessage:", status._message
  print "\tOutput Base URL:", status._baseURL

  if (status._code == 8) or (status._code == 4): # STATUS_DONE || STATUS_FAILED
    break

  print "Waiting 30 seconds"
  sleep(30)
  
  # Query job status
  status = appServicePort.queryStatus(queryStatusRequest(jobID))

# If execution is successful retrieve job outputs. 
if status._code == 8: # 8 = GramJob.STATUS_DONE

  # Create local directory to hold output files
  output_dir = 'meme_out'
  mkdir(output_dir);

  # Instantiate a collection of outputs
  outputs = appServicePort.getOutputs(getOutputsRequest(jobID))

  # Instantiate a request object for retrieving output files
  fileRequest = getOutputAsBase64ByNameRequest()

  # Retreive each output file and save to local directory.
  # MEME writes stdout and stderr to files that aren't listed
  # in the response, so we handle those separately.
  print "\tStandard Output:", outputs._stdOut
  localOutput = open(output_dir + '/' + "stdout.txt", "w")
  fileRequest._jobID = jobID
  fileRequest._fileName = "stdout.txt"
  content = appServicePort.getOutputAsBase64ByName(fileRequest)
  print >>localOutput, content
  localOutput.close()
  print "\tStandard Error:", outputs._stdErr
  localOutput = open(output_dir + '/' + "stderr.txt", "w")
  fileRequest._jobID = jobID
  fileRequest._fileName = "stderr.txt"
  content = appServicePort.getOutputAsBase64ByName(fileRequest)
  print >>localOutput, content
  localOutput.close()
  # Loop for all the other output files.
  for outputFile in outputs._outputFile:
    # List of output files include 'bin' directory which we should skip.
    if outputFile._name == "bin":
      # Bin is a directory, and not of intereest
      # so skip it
      continue
    print "\t" + outputFile._name, ":", outputFile._url
    localOutput = open(output_dir + '/' + outputFile._name, "w")
    fileRequest._jobID = jobID
    fileRequest._fileName = outputFile._name
    content = appServicePort.getOutputAsBase64ByName(fileRequest)
    print >>localOutput, content
    localOutput.close()

