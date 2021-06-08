import numpy as np
import netCDF4 as nc
import pylab as plt
import scipy as sp
import argparse
import os

# command line:
# python plot_SCMjl.py /test/output/BOMEX/data/SolutionRaw/
def main():
    parser = argparse.ArgumentParser(prog='SCM.jl')
    parser.add_argument("folder")
    args = parser.parse_args()
    folder = os.getcwd()+args.folder
    # cwd = os.getcwd()
    # folder = cwd + folder_

    prog_data = nc.Dataset(folder+'/prog_vs_time.nc', 'r')
    aux_data = nc.Dataset(folder+'/aux_vs_time.nc', 'r')
    z = np.array(prog_data.variables['z'])
    time = np.array(prog_data.variables['time'])
    t_end = len(time)
    t_start = np.where((time-3600)>0)[0][0]

    for var in prog_data.variables:
        if not var=='z':
            if not var=='time':
                values = np.array(prog_data.variables[var])
                figname = 'profile_'+var
                plt.figure(figname)
                plt.plot(np.mean(values[t_start:t_end,:], axis = 0), z)
                plt.xlabel(var)
                plt.ylabel('height (m)')
                plt.savefig(folder+'/'+figname+'.pdf')
                plt.close()

                figname = 'contours_'+var
                plt.figure(figname)
                plt.contourf(time, z, np.fliplr(np.rot90(values, k=3)))
                plt.xlabel(var)
                plt.ylabel('height (m)')
                plt.savefig(folder+'/'+figname+'.pdf')
                plt.close()

    for var in aux_data.variables:
        if not var=='z':
            if not var=='time':
                values = np.array(aux_data.variables[var])
                figname = 'profile_'+var
                plt.figure(figname)
                plt.plot(np.mean(values[t_start:t_end,:], axis = 0), z)
                plt.xlabel(var)
                plt.ylabel('height (m)')
                plt.savefig(folder+'/'+figname+'.pdf')
                plt.close()

                figname = 'contours_'+var
                plt.figure(figname)
                plt.contourf(time, z, np.fliplr(np.rot90(values, k=3)))
                plt.xlabel(var)
                plt.ylabel('height (m)')
                plt.savefig(folder+'/'+figname+'.pdf')
                plt.close()

if __name__ == '__main__':
    main()

