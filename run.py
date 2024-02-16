#!/usr/bin/env python3

import argparse
from sys import exit, stderr
from os import path, listdir, mkdir, rename, stat, remove
from timeit import default_timer as timer
from time import strftime, gmtime
from subprocess import call
from shutil import rmtree

def checkEoS():
    if not path.exists(path.abspath("eos")):
        print("Error: no EoS directory. Aborting...", file=stderr)
        exit()
    if not path.exists(path.abspath("eos/eos.py")):
        print("Error: no eos.py file. Aborting...", file=stderr)
        exit()
    call("python3 eos.py --Tmin=0.01 --Tmax=2.00 > eos.dat", shell=True, cwd=path.abspath("eos"))
    if not path.exists(path.abspath("eos/eos.dat")):
        print("Error: no eos file. Aborting...", file=stderr)
        exit()

def extractEvols(collsys, colleng, p, tau0, etas, centrality):
    pParam = f"p{p:.2f}".replace('.', '').replace('-', '')
    if p > 0: pParam += 'p'
    tau0Param = f"tau{tau0:.1f}".replace('.', '')
    evolTarFileDir = path.expanduser("~/TRENTO-VISHNU-EBE")
    evolTarFileDir = path.join(evolTarFileDir, f"TRENTO-VISHNU_{collsys}_{colleng}_etas_{etas}_{pParam}_{tau0Param}")
    evolTarFileDir = path.join(evolTarFileDir, "TProfiles")
    if not path.exists(path.join(evolTarFileDir, f"TProfiles_bin_cent={centrality}.tar.gz")):
        print("Error: no evolutions for fiven parameter set. Aborting...", file=stderr)
        exit()
    call(f"tar -xzf TProfiles_bin_cent={centrality}.tar.gz -C {path.abspath('')}", shell=True, cwd=evolTarFileDir)
    rename(path.abspath(f"TProfiles_bin_cent={centrality}"), path.abspath("tempevols"))

def runPsin(recompile, nThreads, nRange):
    if not path.exists(path.abspath("psin")) or recompile:
        call("g++ cpsin/*.cpp -Wall -fopenmp -O3 -o psin", shell=True, cwd=path.abspath(""))
    eventN = len([f for f in listdir(path.abspath("tempevols")) if "TProfile" in f])
    command  = f"export OMP_NUM_THREADS={nThreads:d}; "
    command += f"./psin --n={nRange} --eventN={eventN:d}"
    call(command, shell=True, cwd=path.abspath(""))

def moveResults(collsys, colleng, p, tau0, etas, centrality):
    resultsDir = path.abspath("results")
    if not path.exists(resultsDir): mkdir(resultsDir)
    resultsDir = path.join(resultsDir, f"results_{collsys}_{colleng}_etas_{etas}_p={p:.2f}_tau0={tau0:.1f}")
    if not path.exists(resultsDir): mkdir(resultsDir)
    resultsDir = path.join(resultsDir, f"results_centrality={centrality}")
    if not path.exists(resultsDir): mkdir(resultsDir)
    fileList = [f for f in listdir(path.abspath("")) if ".dat" in f]
    for aFile in fileList:
        rename(path.abspath(aFile), path.join(resultsDir, aFile))
    if path.exists(path.abspath("tempevols")):
        rmtree(path.abspath("tempevols"))

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

    if args.calculation not in ["psin", "epsn", "avg"]:
        print("Error: calculation type must be on of: psin, epsn, avg. Aborting...", file=stderr)
        exit()

    checkEoS()
    extractEvols(args.collsys, args.colleng, args.p, args.tau0, args.etas, args.centrality)

    if args.calculation == "psin":
        runPsin(args.recompile, args.NUM_THREADS, args.n)
    
    moveResults(args.collsys, args.colleng, args.p, args.tau0, args.etas, args.centrality)