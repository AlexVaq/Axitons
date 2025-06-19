import os, sys
import subprocess
import numpy as np

exe = '../axiton'


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
    exe_source = "../build/axiton"
    run_dir = "test"
    exe_link = os.path.join(run_dir, "axiton")

    os.makedirs(run_dir, exist_ok=True)

    if not os.path.exists(exe_link):
        os.symlink(os.path.abspath(exe_source), exe_link)

    run_dir = "test"
    os.makedirs(run_dir, exist_ok=True)

    # Set initial conditions
    ICgen(1,'axitv',0.3,30, path=run_dir)

    N = 16384
    # Run the binary
    options = ["--size", str(N), "--ctype", "user"]
    command = [exe] + options
    print(command)
    subprocess.run(command, cwd=run_dir)



if __name__ == "__main__":
    run()

