#! /usr/bin/env python

import os
import sys
import optparse

usage = "usage: To be run from DijetRootTreeAnalyzer/ :  python scripts/submit_batch_T2.py -i directory_containing_lists_to_run -o output_directory"

parser = optparse.OptionParser("submitAllGJetsID.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
                  help='run in batch in queue specified as option (default -q cmslong)', 
                  default='cmslong',
                  metavar="QUEUE")

parser.add_option("-i", "--input", dest="input",
                  help="directory containing the lists of samples to be analyzed",
                  )

parser.add_option("-o", "--output", dest="output",
                  help="the directory OUTDIR contains the output of the program",
                  metavar="OUTDIR")

parser.add_option("-m", "--match", dest="match",
                  help="run only the samples containing this string in the name",
                  default="")

parser.add_option('-I', '--interactive',      
                  action='store_true',
                  dest='interactive',      
                  help='run the jobs interactively, 2 jobs at a time',
                  default=False)


(opt, args) = parser.parse_args()
################################################

os.system("mkdir -p batch")
pwd = os.environ['PWD']

if not opt.match:
    os.system("ls "+opt.input+" > config/lists_to_run.txt")
else:
    os.system("ls "+opt.input+" | grep "+opt.match+"  > config/lists_to_run.txt")

ins = open("config/lists_to_run.txt", "r") 
for line in  ins:
    sample = os.path.splitext(line)[0]
    print ("process %s" % sample)
    line = line.rstrip('\n')
    sample = sample.rstrip('\n')

    command = "./main "+opt.input+"/"+line+" config/cutFile_mainDijetSelection.txt dijets/events "+opt.output+"/rootfile_"+sample+" "+opt.output+"/cutEfficiencyFile_"+sample
    print "submit "+command
    print ""
    
    logfile = "logfile_"+sample+".log"
    outputname = "batch/submit_"+sample+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc481\n')
    outputfile.write('cd '+pwd+' \n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write(command+"\n")
    print outputname 
    if opt.interactive==False:
        os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
    else:
        print logfile
        if imc==0: os.system(command+" >&! "+logfile+"&")
        else: os.system(command+" >&! "+logfile)
                
          
