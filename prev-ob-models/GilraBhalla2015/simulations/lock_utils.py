import time
import portalocker # lock files between processes by Jonathan Feinberg

def portalock_open(myfilename):
    myfile = open(myfilename,'r+') # open in rw mode, never write mode
    ## both LOCK_EX and LOCK_SH wait indefinitely for lock to get acquired
    ## i.e. if other process has locked file, wait indefinitely
    ## till lock is released by other process
    portalocker.lock(myfile,portalocker.LOCK_EX) # try to acquire lock
    return myfile

## the problem with writing to a file (mylock below) is that in case:
## one program crashes while files are locked,
## the other program will remain blocked from accessing files!
## can catch exceptions, but portalocker seems fine with one file,
## multiple files seem to be an issue with portalocker.
## go with portalocker and one master lockfile: on process exit, lock ends.
## ensure that you have 'clear\n' in locksimfile.txt if running this process first
def mylock(myfilename,lockstr):
    while True:
        lock_file = portalock_open(myfilename)
        if 'clear' in lock_file.read():
            lock_file.seek(0) # imp: since truncate() removes only from current posn
            lock_file.truncate()
            lock_file.write(lockstr)
            lock_file.flush()
            portalocker.unlock(lock_file)
            lock_file.close()
            break
        portalocker.unlock(lock_file)
        lock_file.close()
        time.sleep(1.0) # don't check too often, wait a second

def myunlock(myfilename):
    lock_file = portalock_open(myfilename)
    lock_file.truncate()
    lock_file.write('clear\n')
    lock_file.flush()
    portalocker.unlock(lock_file)
    lock_file.close()
