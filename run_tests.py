#!/usr/bin/python

# This is the command that needs to be run:
#./build/ECE666/gem5.opt --outdir=<outputDir> ./configs/example/se_splash.py --ruby --num-cpus=<numCPUS> -I <maxInsns> --rootdir=./configs/splash2/codes --benchmark=<benchmark>

from subprocess import Popen, PIPE

#Change this when running regression with a different protocol.  This sets up a separate results dir
PROTOCOL = "MSI"

OUTPUT_DIR = "RunOutputs" + "/" + PROTOCOL
SPLASH_BENCHES = ["Barnes", "OceanContig", "OceanNoncontig", "FMM", "FFT", "Raytrace"]
NUM_CORES = [8, 16, 32, 64]

def buildGem5Command(runDir, numCPUs, maxInsns, benchmark):
  # maxInss == -1 will run the full benchmark
  if maxInsns == -1:
    command = ["./build/ECE666/gem5.opt", "--outdir=" + OUTPUT_DIR + "/" + runDir,
                "./configs/example/se_splash.py", "--ruby", "--num-cpus=" + str(numCPUs), 
                "--rootdir=./configs/splash2/codes", "--benchmark=" + benchmark]
  else:
    command = ["./build/ECE666/gem5.opt", "--outdir=" + OUTPUT_DIR + "/" + runDir,
                "./configs/example/se_splash.py", "--ruby", "--num-cpus=" + str(numCPUs), 
                "-I", str(maxInsns), "--rootdir=./configs/splash2/codes", "--benchmark=" + benchmark]
  return command

def reportError(stdErr):
  errorFound = False
  for line in stdErr.split('\n'):
    if not line.startswith("warn") and (line != ""):
      print line
      errorFound = True
  return errorFound

for numCores in NUM_CORES:
  print "Starting Runs for ", numCores, " Cores..."
  pList = []
  for bench in SPLASH_BENCHES:
    command = buildGem5Command(str(numCores) + "cores_" + bench, numCores, -1, bench)
    print "Starting command: ", command
    pList.append([Popen(command, stdout=PIPE, stderr=PIPE), command])
  for p in pList:
    stdOut, stdErr = p[0].communicate()
    print "Command ", p[1], " finished..."
    reportError(stdErr)

print "Execution Script Finished..."
