import os, sys
import subprocess
import numpy as np


N = 16384
L = 6
nQCD = 7.0

zi = 0.5
zf = 2.0
wDz = 1.5
nstep = 5000


def ICgen(R,ic='flat',sigma=0.1,bst=1,path='test'):
    '''
         the program takes conformal field   theta*R
         and the conformal time derivative   d(theta*R)/d(eta)

         for the moment we just have R=eta  for radiation
                                     R=1    in flat space
    '''
    delta = np.linspace(0,2,1000)
    if ic == 'flat':
        theta = bst*delta/delta
        vheta = theta
        theta = theta*R
    elif ic == 'axitm':
        theta = bst*np.exp(-delta/sigma)
        vheta = theta
        theta = theta*R
    elif ic == 'axitv':
        theta = delta*0
        vheta = bst*np.exp(-delta/sigma)*R
    x = np.column_stack((delta, theta, vheta))
    np.savetxt(path+'/ics.txt', x, delimiter=' ', fmt='%.8e %.8e %.8e')   # X is an array

def run():

    exe_source = "build/src/axiton"
    run_dir = "test"
    exe_link = os.path.join(run_dir, "axiton")

    os.makedirs(run_dir, exist_ok=True)

    if os.path.islink(exe_link):
        if os.readlink(exe_link) != os.path.abspath(exe_source):
            os.remove(exe_link)
            os.symlink(os.path.abspath(exe_source), exe_link)
    else:
        os.symlink(os.path.abspath(exe_source), exe_link)

    # Initial conditions
    ICgen(1,'axitv',0.3,30, path=run_dir)

    # Run the binary
    options = ["--size", str(N), "--ctype", "user", "--zi", str(zi), "--zf", str(zf),
               "--lsize", str(L), "--qcd", str(nQCD), "--wDz", str(wDz), "--steps", str(nstep)]
    command = ["./axiton"] + options
    subprocess.run(command, cwd=run_dir)


if __name__ == "__main__": run()

