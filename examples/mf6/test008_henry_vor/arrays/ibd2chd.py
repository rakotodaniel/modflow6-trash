import numpy as np

#Create a generic function for loading an array
def loadarray(fname, nval=0):
    '''
    Load an array from fname, which can be a file object or a file name.
    if nval = 0, then read all values from the file, otherwise read
    nval values.
    '''
    close_when_done = False
    if isinstance(fname, str):
        f = open(fname, 'r')
        close_when_done = True
    else:
        f = fname
    lst = []
    n = 0
    for line in f:
        ll = line.strip().split()
        for val in ll:
            lst.append(np.float(val))
            n += 1
        if nval > 0:
            if n == nval:
                break
    if close_when_done:
        f.close()
    return np.array(lst)

    
ibound = loadarray('ibound.dat')
strt = loadarray('strt.dat')

f = open('chd.dat', 'w')

 
for n in xrange(ibound.shape[0]):
     if ibound[n] < 0:
         s = '{0} {1} {1}\n'.format(n+1, strt[n])
         f.write(s)
f.close()
