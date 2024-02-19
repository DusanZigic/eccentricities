#!/usr/bin/env python3

import argparse
from sys import exit, stderr
from os import path, listdir, mkdir, rename, remove
from re import findall
import numpy as np
from subprocess import call
from shutil import rmtree, copyfile

def checkEoS():
    if not path.exists(path.abspath("eos")):
        print("Error: no EoS directory. Aborting...", file=stderr)
        exit()
    if not path.exists(path.abspath("eos/eos.py")):
        print("Error: no eos.py file. Aborting...", file=stderr)
        exit()
    if not path.exists(path.abspath("eos/eos.dat")):
        call("python3 eos.py --Tmin=0.01 --Tmax=2.00 > eos.dat", shell=True, cwd=path.abspath("eos"))
    if not path.exists(path.abspath("eos/eos.dat")):
        print("Error: no eos file. Aborting...", file=stderr)
        exit()

def extractEvols(arguments):
    if arguments.calculation in ["plotevol", "eccavgevol"]:
        return
    pParam = f"p{arguments.p:.2f}".replace('.', '').replace('-', '')
    if arguments.p > 0: pParam += 'p'
    tau0Param = f"tau{arguments.tau0:.1f}".replace('.', '')
    evolTarFileDir = path.expanduser("~/TRENTO-VISHNU-EBE")
    evolTarFileDir = path.join(evolTarFileDir, f"TRENTO-VISHNU_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_{pParam}_{tau0Param}")
    evolTarFileDir = path.join(evolTarFileDir, "TProfiles")
    if not path.exists(path.join(evolTarFileDir, f"TProfiles_bin_cent={arguments.centrality}.tar.gz")):
        print("Error: no evolutions for fiven parameter set. Aborting...", file=stderr)
        exit()
    call(f"tar -xzf TProfiles_bin_cent={arguments.centrality}.tar.gz -C {path.abspath('')}", shell=True, cwd=evolTarFileDir)
    if path.exists(path.abspath("tempevols")):
        rmtree("tempevols")
    rename(path.abspath(f"TProfiles_bin_cent={arguments.centrality}"), path.abspath("tempevols"))

def runPsin(arguments):
    if not path.exists(path.abspath("psin")) or arguments.recompile:
        call("g++ cpsin/*.cpp -Wall -fopenmp -O3 -o psin", shell=True, cwd=path.abspath(""))
    eventN = len([f for f in listdir(path.abspath("tempevols")) if "TProfile" in f])
    command  = f"export OMP_NUM_THREADS={arguments.NUM_THREADS:d}; "
    command += f"./psin --n={arguments.n} --eventN={eventN:d}"
    call(command, shell=True, cwd=path.abspath(""))

def runEpsn(arguments):
    psinDir = path.abspath("results")
    psinDir = path.join(psinDir, f"results_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_p={arguments.p:.2f}_tau0={arguments.tau0:.1f}")
    psinDir = path.join(psinDir, f"results_centrality={arguments.centrality}")
    if not path.exists(path.join(psinDir, "psin.dat")):
        print("Error: could not find Psin file. Aborting...", file=stderr)
        exit()
    copyfile(path.join(psinDir, "psin.dat"), path.abspath("psin.dat"))
    if not path.exists(path.abspath("epsn")) or arguments.recompile:
        call("g++ cepsn/*.cpp -Wall -fopenmp -O3 -o epsn", shell=True, cwd=path.abspath(""))
    call(f"export OMP_NUM_THREADS={arguments.NUM_THREADS:d}; ./epsn", shell=True, cwd=path.abspath(""))
    remove(path.abspath("psin.dat"))

def runAvgEvols(arguments):
    psinDir = path.abspath("results")
    psinDir = path.join(psinDir, f"results_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_p={arguments.p:.2f}_tau0={arguments.tau0:.1f}")
    psinDir = path.join(psinDir, f"results_centrality={arguments.centrality}")
    if not path.exists(path.join(psinDir, "psin.dat")):
        print("Error: could not find Psin file. Aborting...", file=stderr)
        exit()
    copyfile(path.join(psinDir, "psin.dat"), path.abspath("psin.dat"))
    if not path.exists(path.abspath("avgevols")) or arguments.recompile:
        call("g++ cavgevols/*.cpp -Wall -fopenmp -O3 -o avgevols", shell=True, cwd=path.abspath(""))
    call(f"./avgevols", shell=True, cwd=path.abspath(""))
    remove(path.abspath("psin.dat"))

def plotEvols(arguments):
    # resultsDir = path.abspath("results")
    # resultsDir = path.join(resultsDir, f"results_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_p={arguments.p:.2f}_tau0={arguments.tau0:.1f}")
    # resultsDir = path.join(resultsDir, f"results_centrality={arguments.centrality}")
    # evolFiles  = [f for f in listdir(resultsDir) if "avgrotevoln" in f]
    # nListFull  = sorted([int(findall(r'\d+', f)[0]) for f in evolFiles])
    # nList      = [n for n in nListFull if n in [1, 2, 3, 4, 5]]
    # for n in nList:
    #     evol = np.loadtxt(path.join(resultsDir, f"avgrotevoln{n:d}.dat"))
    pass

def runEccAvgEvol(arguments):
    if not path.exists(path.abspath("eccavgevols")) or arguments.recompile:
        call("g++ ceccavgevols/*.cpp -Wall -fopenmp -O3 -o eccavgevols", shell=True, cwd=path.abspath(""))
    evolsDir = path.abspath("results")
    evolsDir = path.join(evolsDir, f"results_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_p={arguments.p:.2f}_tau0={arguments.tau0:.1f}")
    evolsDir = path.join(evolsDir, f"results_centrality={arguments.centrality}")
    evolFiles  = [f for f in listdir(evolsDir) if "avgrotevoln" in f]
    for aFile in evolFiles:
        copyfile(path.join(evolsDir, aFile), path.abspath(aFile))
    nList = sorted([int(findall(r'\d+', f)[0]) for f in evolFiles])
    nList = "".join(str(n) for n in nList)
    command  = f"export OMP_NUM_THREADS={arguments.NUM_THREADS:d}; "
    command += f"./eccavgevols --n={nList};"
    call(command, shell=True, cwd=path.abspath(""))
    for aFile in evolFiles:
        remove(path.abspath(aFile))

def moveResults(arguments):
    resultsDir = path.abspath("results")
    if not path.exists(resultsDir): mkdir(resultsDir)
    resultsDir = path.join(resultsDir, f"results_{arguments.collsys}_{arguments.colleng}_etas_{arguments.etas}_p={arguments.p:.2f}_tau0={arguments.tau0:.1f}")
    if not path.exists(resultsDir): mkdir(resultsDir)
    resultsDir = path.join(resultsDir, f"results_centrality={arguments.centrality}")
    if not path.exists(resultsDir): mkdir(resultsDir)
    fileList = [f for f in listdir(path.abspath("")) if ".dat" in f]
    for aFile in fileList:
        rename(path.abspath(aFile), path.join(resultsDir, aFile))
    if path.exists(path.abspath("tempevols")):
        rmtree(path.abspath("tempevols"))
    fileList = [f for f in listdir(path.abspath("")) if ".pdf" in f]
    if len(fileList) > 0:
        figuresDir = path.join(resultsDir, "figures")
        if not path.exists(figuresDir): mkdir(figuresDir)
        for aFile in fileList:
            rename(path.abspath(aFile), path.join(figuresDir, aFile))

if __name__ == "__main__":
    main_dir = path.abspath("")

    parser = argparse.ArgumentParser()
    parser.add_argument('--calculation', type=str,   default = "psin",       help="calculation type: psin, epsn, avg")
    parser.add_argument('--collsys',     type=str,   default = "PbPb",       help="collision system: PbPb or AuAu")
    parser.add_argument('--colleng',     type=str,   default = "5020GeV",    help="collision energy: 5020GeV or 200GeV")
    parser.add_argument('--p',           type=float, default = 0.0,          help="TRENTO p parameter")
    parser.add_argument('--tau0',        type=float, default = 1.0,          help="thermalization time, tau0")
    parser.add_argument('--etas',        type=str,   default = "const",      help="eta/s dependency")
    parser.add_argument('--centrality',  type=str,   default = "30-40%",     help="centrality class")
    parser.add_argument('--n',           type=str,   default = "1-8",        help="range of Fourier expansion coefficients")
    parser.add_argument('--NUM_THREADS', type=int,   default = 50,           help="number of omp threads")
    parser.add_argument('--recompile',   action='store_true', default=False, help="recompile flag")
    args = parser.parse_args()

    if args.calculation not in ["psin", "epsn", "avgevol", "plotevol", "eccavgevol"]:
        print("Error: calculation type must be on of: psin, epsn, avg. Aborting...", file=stderr)
        exit()

    checkEoS()
    extractEvols(args)

    if args.calculation == "psin":
        runPsin(args)
    elif args.calculation == "epsn":
        runEpsn(args)
    elif args.calculation == "avgevol":
        runAvgEvols(args)
    elif args.calculation == "plotevol":
        plotEvols(args)
    elif args.calculation == "eccavgevol":
        runEccAvgEvol(args)
    
    moveResults(args)