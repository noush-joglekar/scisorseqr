#!/usr/bin/env python
# By HT

import os
import pdb
import optparse
import multiprocessing
from multiprocessing import Process, Lock
import sys

################
# global variables
################
lock = Lock();


################
# functions
################
def runCommand(c):
    lock.acquire()
    lock.release()
    os.system(c)
    return


def runcommands(clist, test, pool_size):
    if test == "F":
        if __name__ == '__main__':
            pool = multiprocessing.Pool(processes=pool_size)
            things = pool.map(runCommand, clist)
    else:
        for c in clist:
            print(c)
            print


################
# main
################
def main():
    print("## 0. params")
    print("# 0a. read params")
    parser = optparse.OptionParser()
    parser.add_option("--commandFile", action="store", type="string", dest="commandFile", help="a file to read")
    parser.add_option("--n", action="store", type="int", dest="n", help="a file to read")

    (option, args) = parser.parse_args()

    print("# 0b. print params")
    print("commandFile=" + option.commandFile)
    print("n=" + str(option.n))

    print("## 1. reading organizational data")
    commandList = []
    f = open(option.commandFile, 'r')
    for line in f:
        commandList.append(line.rstrip())

    print("## 2. execution")
    runcommands(commandList, "F", option.n)


################
# execution
################
main()
