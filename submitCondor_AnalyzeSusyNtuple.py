#!/usr/bin/env python

import sys
import optparse
import commands
import os
import tarfile
import glob
import time

#######################
# Get options
#######################

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")
parser.add_option ('--prefix', dest='prefix', type='string',
                   default = 'NONE',
                   help="directory containing input susy ntuples")
parser.add_option ('--o', dest='baseOutdir', type='string',
                   default = 'dummy',
                   help="base out directory")
parser.add_option ('--jdl', dest='jdlBase', type='string',
                   default = 'condor_TEMPLATE.jdl',
                   help="condor jdl base")
parser.add_option ('--filelist_input',dest='filelist_input', type='string',
                   default = 'NONE')
parser.add_option ('--script', dest='scriptBase', type='string',
		   default = 'condor_BASH.sh',
		   help="condor executable script base")
parser.add_option ('--ana', dest='anaBase', type='string',
		   default = 'ana.C',
		   help="root macro, your ana.C")
parser.add_option ('--njobs', dest='njobs', type='int',
                   default = '-1',
                   help="number of jobs, default = 1 job per file")
parser.add_option ('--test', action="store_true",
                   dest="test", default=False,
                   help="Just testing")
parser.add_option ('--json', dest='json', type='string',
                   default = '',
                   help="JSON file used to select events")
parser.add_option ('--inputFolders', dest='inputFolders', type='string',
		   default='',
		   help='Comma-separated list of folders containing nTuple files')

options, args = parser.parse_args()

prefix = options.prefix
test = options.test
baseOutdir = options.baseOutdir
jdlBase = options.jdlBase
filelist_input = options.filelist_input
scriptBase = options.scriptBase
anaBase = options.anaBase
njobs = options.njobs
json = options.json
inputFolders = options.inputFolders

cwd = os.getcwd()

if not os.path.isdir(baseOutdir) :
    os.system("mkdir -p "+baseOutdir)
    print "Making directory %s." % baseOutdir
else:
    print "Output directory %s already exists.  Exiting." % baseOutdir
    sys.exit()

FILE_LIST = open(filelist_input,'r')

allFiles_tmp = []
print inputFolders
print "Prefix: %s" %prefix
if len(inputFolders) != 0:
    for folder in inputFolders.split(':'):
        if folder[-1] != '/':
            folder = folder+'/'
	if filelist_input !='NONE':
            for line in FILE_LIST:
                print line.strip()[-4:]
                if line.strip()[-4:]=="root":
                    line_noR=line.rstrip('\n')
                    allFiles_tmp.append(folder+line_noR)
        else:
            allFiles_tmp += glob.glob(folder+"*.root")
        
else:
    print "Zero files given, already done! Exiting."
    sys.exit()
allFiles = []    
if prefix !="NONE":
    for files in allFiles_tmp:
        allFiles.append(prefix+files)
        print files
else:
    allFiles=allFiles_tmp

nfiles = len(allFiles)
print "Number of files: %d" %nfiles
# If -1, make 1 job per file:
if njobs == -1: njobs = nfiles

# Split files into jobs as evenly as possible:

baseFilesPerJob = nfiles/njobs
extraFiles = nfiles%njobs

numberOfFiles = {}
for ijob in range(njobs):
    numberOfFiles[ijob] = baseFilesPerJob

for ijob in range(extraFiles):
    numberOfFiles[ijob] += 1


files_for_jobs = {}
IFILE = 0
for ijob in range(njobs):
    nfiles_ijob = numberOfFiles[ijob]
    files_for_jobs[ijob] = []
    for ifile in range(nfiles_ijob):
        files_for_jobs[ijob].append( allFiles [IFILE] )
        IFILE += 1

if IFILE != nfiles:
    print "Mismatch.", IFILE, nfiles-1,"  Exiting."
    sys.exit()

for ijob in range(njobs):

    jobid = str(ijob)

    filelistname = "filelist_"+jobid
    filelist = open(baseOutdir+"/"+filelistname, 'w')

    for ifile in files_for_jobs[ijob]:
        filelist.write(ifile+"\n")
    filelist.close()
if not os.path.isdir(baseOutdir+"/JobOut"): os.system("mkdir "+baseOutdir+"/JobOut")
#print "Working direcotry for filelist"
#os.system("pwd")
#os.system("ls -l -h testJob/filelist_0")
path_jdl = baseOutdir+"/"+jdlBase
path_script = baseOutdir+"/"+scriptBase
path_ana = baseOutdir+"/"+anaBase

commandList = []
commandList.append("cp "+jdlBase+" "+baseOutdir)
commandList.append("cp "+json+" "+baseOutdir)
#commandList.append("cp "+filelistname+" "+baseOutdir)
commandList.append("cp "+scriptBase+" "+baseOutdir)
commandList.append("cp "+anaBase+" "+baseOutdir)

commandList.append('replace NJOBS  '+str(njobs)+' -- '+path_jdl)
commandList.append('replace SCRIPT  '+scriptBase+' -- '+path_jdl)
commandList.append('replace ANALYZER  '+anaBase+' -- '+path_jdl)
#commandList.append("replace FILELIST '+filelistname+' -- "+path_jdl)
commandList.append('replace ANALYZER  '+anaBase+' -- '+path_script)

if json == '':
    commandList.append('replace JSON  " " -- '+path_jdl)
    commandList.append('replace JSON  " " -- '+path_script)
    commandList.append('replace JSON  " " -- '+path_ana)
else:
    commandList.append('replace JSON  '+json+' -- '+path_jdl)
    commandList.append('replace JSON  '+json+' -- '+path_script)
    commandList.append('replace JSON  '+json+' -- '+path_ana)

for command in commandList:
    os.system(command)
if not test:
#    os.system("pwd")
#    os.system("ls -l -h filelist_0")
    os.chdir(baseOutdir)
#    os.system("pwd")
#    os.system("ls")
#    os.system("ls -l -h filelist_0")
    os.system("tar -czf fileLists.tgz filelist_*")
    os.system("condor_submit "+jdlBase)
    os.chdir(cwd)





