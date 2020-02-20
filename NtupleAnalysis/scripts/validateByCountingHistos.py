#!/usr/bin/env python

import os
import sys
import re

import ROOT

root_re = re.compile("(?P<rootfile>([^/]*))\.root")

sAllHistos = set()
sEmptyHistos = set()

def execute(cmd):
    f = os.popen4(cmd)[1]
    ret=[]
    for line in f:
        ret.append(line.replace("\n", ""))
    f.close()
    return ret

def ls(path=""):
    global sAllHistos
    global sEmptyHistos

    ROOT.gFile.cd(path)
    gDir = ROOT.gFile.CurrentDirectory()
    print "dir",path
    keys = gDir.GetListOfKeys()
    for i in range(gDir.GetNkeys()):
        keyname = keys.At(i).GetName()
        print "check key",keyname
        obj = gDir.Get(keyname)

        if isinstance(obj, ROOT.TH1F) or isinstance(obj, ROOT.TH2F):
            sAllHistos.add(os.path.join(path,keyname))
            if not lsHisto(obj):
                sEmptyHistos.add(os.path.join(path,keyname))

        #if isinstance(obj, ROOT.TProfile):
        #    lsHisto(obj)

        #if isinstance(obj, ROOT.TTree):
        #    lsTree(obj)

        if isinstance(obj, ROOT.TDirectoryFile):
            ls(os.path.join(path,keyname))

    ROOT.gFile.cd(path)


def lsHisto(histo):
    gDir = ROOT.gFile.CurrentDirectory()
    dirName = gDir.GetName()
    if histo.GetEntries() > 0:
        return True
    return False

def lsTree(tree):
    print "Tree",tree.GetName(),", entries: ",tree.GetEntries(),", branches:"

    for branch in tree.GetListOfBranches():
        print "    ",branch.GetName()

returnFiles = []
def findFiles(dirs):
    global returnFiles
    for d in dirs:
        cands = execute("ls %s"%d)
        for c in cands:
            filepath = c
            if os.path.isdir(d):
                filepath = os.path.join(d,c)
            if os.path.isdir(filepath):
                findFiles([filepath])
            else:
                returnFiles.append(filepath)
    return returnFiles

def main():

    if len(sys.argv) == 1:
        print "\n"
        print "### Usage:   ",os.path.basename(sys.argv[0])," <dir or root file>\n"
        print "\n"
        sys.exit()

    files = findFiles(sys.argv[1:])
    for file in files:

        i = 2
        tree  = ""
        while i < len(sys.argv):
            if i != len(sys.argv)-1 and (sys.argv[i] == "-t" or sys.argv[i] == "-tree"):
                tree = sys.argv[i+1] 
            i = i + 1

        match = root_re.search(file)

        sAllHistos.clear()
        sEmptyHistos.clear()
        if match and os.path.isfile(file):
            fIN = ROOT.TFile.Open(file,'r')

            ls()

            print "File",file
            print "    Number of all histograms   ",len(sAllHistos)
            print "    Number of empty histograms ",len(sEmptyHistos)
            if len(files) == 1:
                print "    Empty histograms:"
                for h in sEmptyHistos:
                    print "       ",h
                print "    Number of all histograms   ",len(sAllHistos)
                print "    Number of empty histograms ",len(sEmptyHistos)

if __name__ == "__main__":
    main()
