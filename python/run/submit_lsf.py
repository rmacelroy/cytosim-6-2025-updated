#!/usr/bin/env python3
#
# A script to submit jobs to the Platform LSF queuing system
#
# F. Nedelec, 10.2007 -- 12.2016
# Adapted for LSF by Andre Claude Clapson, November 2012

"""
    Submit an array of jobs to the Platform LSF queue, to be started by 'go_sim.py'
    
Syntax:
    
    submit.py ARG [mem=????] [avx=1] [queue=????] file1 [file2] [file3] [...]
    
    ARG is passed as arguments to 'go_sim.py' (please check go_sim.py help),
    and you must use quotes if you have multiple arguments to group them together.
    
    Unless specified otherwise, the queue is 'medium_priority'.
    The amount of requested memory (default=2G) can be specified as:
       mem=1024 (for 1 GB)
       mem=512  (for 512 MB)
       ...
    
Examples:
    
    1. submit.py sim config.cym
       Submit one job to run the config file provided
    
    2. submit.py sim config1.cym config2.cym config3.cym ...
       Submit one job for each config file provided. They will run in parallel.

    3. submit.py 'sim preconfig.py' config.cym.tpl
       Submit one job where the config file will be generated by preconfig.py from config.cym.tpl
       The executable will run sequentially on one node, for each file that is generated.
    
    4. submit.py 'sim 10' config.cym
       Submit one job to run the same config file 10 times sequentially on one node.
    
    5. submit.py 'sim 10' config1.cym config2.cym config3.cym ...
       Submit a different job for each config file provided, with 10 repeats of each
    
    
Last updated 30 July 2013
F. Nedelec
"""


import sys, os, shutil, subprocess


out    = sys.stderr
bsub   = 'bsub'
queue  = 'medium_priority'
mem    = '3096'
jdir   = 'job00'

#------------------------------------------------------------------------

def makeDirectory(name):
    """create a directory with specified name, if it does not exists"""
    if os.path.isdir(name):
        if not os.listdir(name):
            return 1
    else:
        try:
            os.mkdir(name)
            return 2
        except OSError:
            out.write("Error: directory '%s' could not be created\n" % path)
            sys.exit()
    return 0


def makeEmptyDirectories(names):
    """Create a set of empty directories that start with the given names"""
    cnt = 0
    done = []
    while cnt < 10000:
        add = str(cnt)
        res = []
        for p in names:
            res.append(makeDirectory(p+add))
        done = sum(map(lambda x:x>0, res))
        if done == len(names):
            return add
        #remove directories that were made:
        for p in names:
            if res.pop(0) == 2:
                os.rmdir(p+add)
        cnt += 1
    return ''


def makeNumberedDirectory(root):
    """Create an empty directories that start with 'root'"""
    cnt = 0
    while cnt < 10000:
        name = root + '%02i' % cnt
        try:
            os.mkdir(name)
            return name
        except OSError:
            cnt += 1
    out.write("Error: directory '%s' could not be created\n" % root)
    return ''

#------------------------------------------------------------------------

def makeScript(filename, cmd):
    """create an executable file containing the commands"""
    file = open(filename, 'w')
    file.write(cmd)
    file.close()
    os.chmod(filename, 0700)


def job(cwd, conf, jarg):
    """in a string, return the commands necessary to run go_sim.py"""
    cmd  = '#!/bin/bash\n'
    # cmd += 'set -x;\n'
    cmd += 'cd %s;\n' % cwd
    cmd += 'touch %s\n' % conf
    # the job will call go_sim.py once:
    cmd += './go_sim.py %s %s;\n' % (jarg, conf)
    cmd += 'mv '+conf+' '+jdir+'/done/.;\n'
    cmd += '\\rm $0;\n'
    # cmd += 'mv $0 '+jdir+'/done/.;\n'
    return cmd


def sub(exe):
    """return a command that will submit one job"""
    # specify memory, shell, minimum number of cores and queue
    cmd  = [bsub, '-L', '/bin/bash', '-n', '1', '-q', queue]
    cmd += ['-M', mem, '-R', 'rusage[mem=%s]' % mem]
    # redirect stderr and sdtout to files:
    cmd += ['-oo', jdir+'/logs/out_%J.txt']
    cmd += ['-eo', jdir+'/logs/err_%J.txt']
    # executable
    cmd += [exe]
    return cmd


def subArray(first, last, todo):
    """return a command that will submit a job-array"""
    # set memory , queue , CPU number and request logon on host
    cmd  = [bsub, '-L', '/bin/bash', '-n', '1', '-q', queue]
    cmd += ['-M', mem, '-R', 'rusage[mem=%s]' % mem]
    # defines the jobarray name and index range
    cmd += ['-J', '%s[%i-%i]' % (todo, first, last)]
    # redirect stderr and sdtout to files:
    cmd += ['-oo', jdir+'/logs/out_%I.txt']
    cmd += ['-eo', jdir+'/logs/err_%I.txt']
    # job script
    cmd += [todo+'$LSB_JOBINDEX']
    return cmd

#------------------------------------------------------------------------

def execute(cmd):
    """execute given command with subprocess.call()"""
    try:
        val = subprocess.call(cmd)
        if val:
            out.write("Error: command failed with value %i\n" % val)
            print(cmd)
    except OSError:
        out.write("Error: command `%s' failed:\n" % cmd[0])
        print(cmd)


def main(args):
    """submit jobs, depending on the arguments provided"""
    global bsub, mem, queue, jdir
    
    #find bsub command:
    proc = subprocess.Popen(['which', 'bsub'], stdout=subprocess.PIPE)
    if proc.wait():
        out.write("Error: LSF command `bsub' not found\n")
    #sys.exit()
    else:
        bsub = proc.stdout.readline().strip()
        #print('|'+bsub+'|')
    
    # first argument is used for go_sim.py:
    jarg = args.pop(0)
    
    # this is to catch old-style invocation
    if jarg.endswith(".cym"):
        out.write("Error: you must first provide arguments to go_sim.py\n")
        out.write("       enter `submit.py help` for info\n")
        sys.exit()
    
    # make new directories for this job
    jdir = makeNumberedDirectory('job')
    os.mkdir(os.path.join(jdir, 'todo'))
    os.mkdir(os.path.join(jdir, 'done'))
    os.mkdir(os.path.join(jdir, 'save'))
    os.mkdir(os.path.join(jdir, 'logs'))

    print("    go_sim.py will run `%s' in %s" % (jarg, jdir))

    avx  = 0
    cnt  = 1
    id   = 0
    file = ''
    cwd  = os.getcwd()
    todo = os.path.join(jdir,'todo')
    
    for arg in args:
        if os.path.isfile(arg):
            if os.access(arg, os.X_OK):
                out.write("Error: file `%s' should not be executable\n" % arg)
                sys.exit()
            for x in range(cnt):
                id += 1
                conf = todo + '/config%04i.cym' % id
                shutil.copyfile(arg, conf)
                file  = todo + '/R' + str(id)
                jargs = jarg + ' name=run%04i'%id + ' park='+jdir+'/save'
                makeScript(file, job(cwd, conf, jargs))
        elif arg.isdigit():
            cnt = int(arg)
        elif arg.startswith('mem='):
            mem = arg[4:]
        elif arg.startswith('avx='):
            avx = int(arg[4:])
        elif arg.startswith('queue='):
            queue = arg[6:]
        else:
            out.write("Error: I do not understand argument `%s'\n" % arg)
            sys.exit()

    if id == 0:
        out.write("Error: you need to specify at least one config file\n" % arg)
        sys.exit()
    
    if id > 1:
        print("    submit.py created %i scripts in `%s'" % (id, todo))
        cmd = subArray(1, id, todo+'/R')
    else:
        cmd = sub(file)
    # request nodes with support for Intel AVX technology:
    if avx:
        cmd += ['-m', 'intelavx']
    execute(cmd)


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1].endswith("help"):
        print(__doc__)
    else:
        main(sys.argv[1:])

