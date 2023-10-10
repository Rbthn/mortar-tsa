import numpy as np
import matplotlib.pyplot as plt
import os
import gmsh
import csv
from mpl_toolkits.axes_grid1 import host_subplot

main_dir = os.path.dirname(__file__)
os.chdir(main_dir)


def save_csv(file, labels, data):
    with open(file, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(labels)
        for r in data:
            writer.writerow(r)
    return

def reset_matheval():
    gmsh.plugin.setNumber('MathEval', 'View', -1)
    gmsh.plugin.setNumber('MathEval', 'OtherView', -1)
    gmsh.plugin.setNumber('MathEval', 'TimeStep', -1)
    gmsh.plugin.setNumber('MathEval', 'OtherTimeStep', -1)
    gmsh.plugin.setNumber('MathEval', 'ForceInterpolation', 0)
    gmsh.plugin.setNumber('MathEval', 'PhysicalRegion', -1)
    gmsh.plugin.setString('MathEval', 'Expression0', '')
    gmsh.plugin.setString('MathEval', 'Expression1', '')
    gmsh.plugin.setString('MathEval', 'Expression2', '')
    gmsh.plugin.setString('MathEval', 'Expression3', '')
    gmsh.plugin.setString('MathEval', 'Expression4', '')
    gmsh.plugin.setString('MathEval', 'Expression5', '')
    gmsh.plugin.setString('MathEval', 'Expression6', '')
    gmsh.plugin.setString('MathEval', 'Expression7', '')
    gmsh.plugin.setString('MathEval', 'Expression8', '')

'''
Start gmsh and logger
'''
def start_gmsh():
    gmsh.initialize(['-v', '5', '-setnumber', 'General.Terminal', '0'])
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()
    return


'''
Write logs, stop gmsh
'''
def stop_gmsh(logfile):
    gmsh.logger.stop()
    gmsh.finalize()
    return


'''
Compute L2-Norm of local euclidian error
'''
def compute_err_L2(quantity_file_name, reference_file_name, logfile):
    # Hide initial "Increasing process stack size" message
    start_gmsh()

    has_reference = reference_file_name != ""

    if has_reference:
        gmsh.open(reference_file_name) # read analytical / reference solution. View 0
    gmsh.open(quantity_file_name) # read mortar/TSA solution. View 1

    # Use MathEval to compute sq. euclidian error
    tags = gmsh.view.getTags()
    idx = [gmsh.view.getIndex(t) for t in tags]
    gmsh.plugin.setNumber('MathEval', 'View', idx[0]) # use reference solution as main view ('vi')

    if has_reference:
        gmsh.plugin.setNumber('MathEval', 'OtherView', idx[1]) # use solution as other view ('wi')
        gmsh.plugin.setString('MathEval', 'Expression0', '(v0-w0)^2+(v1-w1)^2+(v2-w2)^2') # Squared error. View 2
    else:
        gmsh.plugin.setString('MathEval', 'Expression0', '(v0)^2+(v1)^2+(v2)^2')
    gmsh.plugin.run('MathEval')
    reset_matheval()

    # Use Integrate to compute sq. L2 error from local euclidian error
    gmsh.plugin.run('Integrate')
    file_name = "tmp.csv"
    gmsh.view.write(gmsh.view.getTags()[-1], file_name)

    # read sq. L2 error results, compute sqrt.
    # Doing this step in python to get a more convenient output
    sq_L2_data = np.genfromtxt(file_name, delimiter=' ', ndmin=2)
    # i, ti, ?, ?, x, y, z, L2
    timesteps = sq_L2_data[:, 1]
    sq_L2 = sq_L2_data[:, 7]
    L2 = np.sqrt(sq_L2)

    stop_gmsh(logfile)
    return timesteps, L2


I = np.genfromtxt("I.csv", delimiter=' ', ndmin=2)
I_ref = np.genfromtxt("I_ref.csv", delimiter=' ', ndmin=2)
V = np.genfromtxt("V.csv", delimiter=' ', ndmin=2)
V_ref = np.genfromtxt("V_ref.csv", delimiter=' ', ndmin=2)

###
### plot voltages over time
###

host = host_subplot(111)
plt.xlim(0, 7)
ax_i = host.twinx()

host.plot(V_ref[:,0], -1000*V_ref[:,1], color='#440154', label='Reference solution voltage')
host.plot(V[:,0], 1000*V[:,1], color='#21918c', ls='--', label='Mortar-TSA voltage')
ax_i.plot(I_ref[:, 0], I_ref[:, 1], color='#31688e', linestyle='-.', label='Source current')

host.set_xlabel("Time [s]")
host.set_ylabel("Voltage [mV]")
ax_i.set_ylabel("Current [A]")
host.legend(labelcolor="linecolor")
host.grid(True)
host.set_ylim(-0.1, 5.5)
ax_i.set_ylim(-1, 55)

# plt.savefig("h-phi-voltage.pdf")
plt.show()

### plot L2 error of B over time
'''
time, ref_L2 = compute_err_L2("b_ref.pos", "", "")
time_abs, abs_L2 = compute_err_L2("b.pos", "b_ref.pos", "")
rel_L2 = abs_L2 / ref_L2

labels = ["time", "reference L2", "absolute L2", "relative L2"]
data = np.stack((time, ref_L2, abs_L2, rel_L2), axis=1)
save_csv("test_case_error.csv", labels, data)
'''

L2_data = np.genfromtxt("test_case_error.csv", delimiter=',', ndmin=2, skip_header=1)
time = L2_data[:, 0]
ref_L2 = L2_data[:, 1]
abs_L2 = L2_data[:, 2]
rel_L2 = L2_data[:, 3]

fig_b, ax_b = plt.subplots()

ax_b.semilogy(time, ref_L2, color='#440154', label='Reference solution')
ax_b.semilogy(time, abs_L2, color='#21918c', linestyle='--', label='Absolute error')
ax_b.semilogy(time, rel_L2, color='#31688e', linestyle='-.', label='Relative error')
plt.xlim(0, 7)
plt.ylim(1e-13, 1)
ax_b.yaxis.get_major_locator().set_params(numticks=999)
ax_b.yaxis.get_minor_locator().set_params(numticks=999, subs=[2, 4, 6, 8])

ax_b.set_xlabel("Time [s]")
ax_b.set_ylabel("$L^2$-Norm")
ax_b.legend(loc='best')
ax_b.grid(True, which='both', lw=0.5)
#plt.savefig("h-phi-b-err.pdf")
plt.show()

print('done!')