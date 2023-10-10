#!/bin/python3
#
import sys, os
import subprocess, resource, datetime # time and log subprocess
import platform, psutil # system info
import gmsh, csv, numpy as np # error calculation
import argparse # process CLI options for this script
from glob import glob # use shell wildcards when saving output files

###
### Helper functions
###

def log_str_pre(command: str, isodate: str):
    return "#"*80+"\n"+"#"*80+"\n### Running command:\t%s\n### at\t\t\t%s\n###\n" % (command, isodate)
def log_str_post(isodate: str, elapsed: float):
    return "\n###\n### Finished at\t\t%s\n### CPU time:\t\t%f s\n" % (isodate, elapsed) +"#"*80+"\n"+"#"*80+"\n"*3
def isodate():
    return datetime.datetime.now().strftime("%Y-%m-%dT%H_%M_%S%z")
# format: display name, reference L2, absolute error, relative error
def iteration_report(iteration, labels, data):
    ret = "###\tIteration\t"+str(iteration)
    if (len(data) != 0):
        ret += "\n###\tName"
        for label in labels:
            ret += "\t" + label
    for line in data:
        ret += "\n###"
        for i in range(len(line)):
            if i == 0:
                ret += "\t" + line[i]
            else:
                ret += "\t%f" % (line[i])
    ret += '\n'
    return ret

def time_and_log_subprocess(command: str, args: list, stdout_filename: str, stderr_filename: str):
    # Print info pre
    start_isodate = isodate()
    command_str = command
    for ele in args:
        command_str = command_str + " " + ele

    print_and_log(log_str_pre(command_str, start_isodate), 1, stdout_filename)

    ###
    ### actual subprocess.
    ### If stdout/err files cannot be opened, output to terminal
    ###
    args.insert(0, command) # put command at the front of arg list
    end_isodate = ""
    cpu_time = -1
    return_code = -1


    try:
        stdout_f = open(stdout_filename, 'a')
        stderr_f = open(stderr_filename, 'a')
    except IOError:
        start = resource.getrusage(resource.RUSAGE_CHILDREN)
        process = subprocess.run(args)
        return_code = process.returncode
        end = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_time = end.ru_utime - start.ru_utime
    else:
        with stdout_f, stderr_f:
            start = resource.getrusage(resource.RUSAGE_CHILDREN)
            process = subprocess.run(args, stdout=stdout_f, stderr=stderr_f)
            return_code = process.returncode
            end = resource.getrusage(resource.RUSAGE_CHILDREN)
            cpu_time = end.ru_utime - start.ru_utime
            end_isodate = isodate()

    # Print info post
    print_and_log(log_str_post(end_isodate, cpu_time), 1, stdout_filename)
    assert return_code == 0
    return cpu_time

'''
Collecting system information, such as:
- CPU
- Mem total
- Mem free
- Hostname
'''
def collect_system_info(stdout_filename=''):
    def get_size(bytes, suffix="B"):
        """
        Scale bytes to its proper format
        e.g:
            1253656 => '1.20MB'
            1253656678 => '1.17GB'
        """
        factor = 1024
        for unit in ["", "K", "M", "G", "T", "P"]:
            if bytes < factor:
                return f"{bytes:.2f}{unit}{suffix}"
            bytes /= factor


    uname = platform.uname()
    cpu_freq = psutil.cpu_freq()
    swap = psutil.swap_memory()
    vmem = psutil.virtual_memory()
    info_str = """###
### System
###

System:\t\t\t%s
Hostname:\t\t%s
Architecture:\t\t%s
Release:\t\t%s


###
### CPU
###

Physical cores:\t\t%s
Logical cores:\t\t%s
CPU Freq. Min:\t\t%.2f MHz
CPU Freq. Max:\t\t%.2f MHz
CPU Freq. Current:\t%.2f MHz


###
### MEM
###

MEM total:\t\t%s
MEM used:\t\t%s
MEM available:\t\t%s
Swap total:\t\t%s
Swap used:\t\t%s\n\n""" % (
        uname.system,
        uname.node,
        uname.machine,
        uname.release,
        psutil.cpu_count(logical=False),
        psutil.cpu_count(logical=True),
        cpu_freq.min,
        cpu_freq.max,
        cpu_freq.current,
        get_size(vmem.total),
        get_size(vmem.used),
        get_size(vmem.available),
        get_size(swap.total),
        get_size(swap.used)
    )

    return print_and_log(info_str, 1, stdout_filename)


'''
Write text to logfile. If current VERBOSITY is greater or equal to verbosity_threshold, print to terminal as well.
'''
def print_and_log(text, verbosity_threshold, logfile):
    if VERBOSITY >= verbosity_threshold:
        print(text)
    with open(logfile, 'a') as file:
        if text is list:
            for line in text:
                file.write(text)
        else:
            file.write(text)
        return 0
    return 1



###
### Error calculation
###

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
    return L2

'''
Compute L1 Norm of local euclidian error
'''
def compute_err_L1(quantity_file_name, reference_file_name, logfile):
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
        gmsh.plugin.setString('MathEval', 'Expression0', 'Sqrt((v0-w0)^2+(v1-w1)^2+(v2-w2)^2)') # Absolute error. View 2
    else:
        gmsh.plugin.setString('MathEval', 'Expression0', 'Sqrt(v0^2+v1^2+v2^2)') # Norm. View 1
    gmsh.plugin.run('MathEval')
    reset_matheval()

    # Use Integrate to compute L1 error from local euclidian error
    gmsh.plugin.run('Integrate')
    file_name = "tmp.csv"
    gmsh.view.write(gmsh.view.getTags()[-1], file_name)

    # read L1 error results
    L1_data = np.genfromtxt(file_name, delimiter=' ', ndmin=2)
    # i, ti, ?, ?, x, y, z, L2
    timesteps = L1_data[:, 1]
    L1 = L1_data[:, 7]

    stop_gmsh(logfile)
    return L1

'''
Compute LInf (=Max) Norm of local euclidian error
'''
def compute_err_LInf(quantity_file_name, reference_file_name, logfile):
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
        gmsh.plugin.setString('MathEval', 'Expression0', 'Sqrt((v0-w0)^2+(v1-w1)^2+(v2-w2)^2)') # Absolute error. View 2
    else:
        gmsh.plugin.setString('MathEval', 'Expression0', 'Sqrt((v0)^2+(v1)^2+(v2)^2)')
    gmsh.plugin.run('MathEval')
    reset_matheval()

    # Use MixMax to compute maximum error
    gmsh.plugin.setNumber('MinMax', 'Argument', 1)
    gmsh.plugin.run('MinMax')
    max_file_name = "tmp.csv"
    gmsh.view.write(gmsh.view.getTags()[-1], max_file_name)

    # read L1 error results
    LInf_data = np.genfromtxt(max_file_name, delimiter=' ', ndmin=2)
    # i, ti, ?, ?, x, y, z, Max
    timesteps = LInf_data[:, 1]
    LInf = LInf_data[:, 7]

    stop_gmsh(logfile)
    return LInf


'''
Save error results in csv
'''
def save_csv(file, labels, data):
    with open(file, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(labels)
        for r in data:
            writer.writerow(r)
    return


'''
Start gmsh and logger
'''
def start_gmsh():
    gmsh.initialize(['-v', '5', '-setnumber', 'General.Terminal', '0'])
    if VERBOSITY < 2:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.logger.start()
    return


'''
Write logs, stop gmsh
'''
def stop_gmsh(logfile):
    logs = gmsh.logger.get()
    with open(logfile, 'a') as file:
        for line in logs:
            file.write(line + '\n')
    gmsh.logger.stop()
    gmsh.finalize()
    return


###
### MAIN
###

# Termination criteria, settings
TOL_REL_ERR = 5e-3
TOL_ABS_ERR = 1e-8
MAX_IT = 5
REFERENCE_REFINEMENT = 2
LOG_NAME = "test.log"
VERBOSITY = 0
QUANTITIES = []

parser = argparse.ArgumentParser(
    prog="Main script",
    description="",
)

parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-q', '--quantity', nargs='+', type=str)
parser.add_argument('--tol_rel', type=float)
parser.add_argument('--tol_abs', type=float)
parser.add_argument('--max_it', type=int)
parser.add_argument('--logfile', type=str)
parser.add_argument('--reference_refine', type=int)
parser.add_argument('--analytical', action='store_true')
parser.add_argument('--steady_state', action='store_true')
parser.add_argument('target_dir_rel')

args = parser.parse_args()
VERBOSITY = args.verbose
USE_ANALYTICAL = args.analytical
STEADY_STATE = args.steady_state
if args.quantity:
    QUANTITIES = args.quantity
if args.tol_rel:
    TOL_REL_ERR = args.tol_rel
if args.tol_abs:
    TOL_ABS_ERR = args.tol_abs
if args.max_it:
    MAX_IT = args.max_it
if args.logfile:
    LOG_NAME = args.logfile
if args.reference_refine:
    REFERENCE_REFINEMENT = args.reference_refine


target_dir_rel = args.target_dir_rel # problem dir, e.g. code/h-phi_3d. Relative from script call location.
main_dir = os.path.dirname(__file__) # main dir, i.e. code/common

# navigate to target
os.chdir(target_dir_rel)
target_dir_abs = os.getcwd()
logfile = os.path.join(target_dir_abs, LOG_NAME)

# info
collect_system_info(logfile)

# if not specified, determine reference refinement level
if not USE_ANALYTICAL:
    if not args.reference_refine:
        print_and_log("No reference refinement level supplied.", 0, logfile)
        print_and_log("Trying to determine a reasonable level...", 0, logfile)
        # solve with refine 0
        time_mesh_ref    = time_and_log_subprocess(command='python3',  args=['mesh.py','--reference', '--refine', '0'], stdout_filename=logfile, stderr_filename=logfile)
        solve_args = ['solve.py','--reference']
        if STEADY_STATE:
            solve_args.append('--steady_state')
        time_solve_ref   = time_and_log_subprocess(command='python3', args=solve_args, stdout_filename=logfile, stderr_filename=logfile)

        REFERENCE_REFINEMENT = 1
        while True:
            # rename old results
            for q in QUANTITIES:
                time_and_log_subprocess(command='mv', args=[q+"_ref.pos", q+"_ref_old.pos"], stdout_filename=logfile, stderr_filename=logfile)
            # solve for current refine
            time_mesh_ref    = time_and_log_subprocess(command='python3',  args=['mesh.py','--reference', '--refine', str(REFERENCE_REFINEMENT)], stdout_filename=logfile, stderr_filename=logfile)
            solve_args = ['solve.py','--reference']
            if STEADY_STATE:
                solve_args.append('--steady_state')
            time_solve_ref   = time_and_log_subprocess(command='python3', args=solve_args, stdout_filename=logfile, stderr_filename=logfile)

            # determine error. Use current refine as reference.
            quantities_ref_L2 = [compute_err_L2(q+"_ref.pos", "", logfile)[-1] for q in QUANTITIES]
            quantities_abs_L2 = [compute_err_L2(q+"_ref_old.pos", q+"_ref.pos", logfile)[-1] for q in QUANTITIES]
            quantities_rel_L2 = [quantities_abs_L2[i] / quantities_ref_L2[i] for i in range(len(QUANTITIES))]

            logtext = "\tReference refinement level="+str(REFERENCE_REFINEMENT)+": relative L2 error " + str(quantities_rel_L2)
            print_and_log(logtext, 0, logfile)

            rel_tol_fulfilled = list(filter(lambda x: x < TOL_REL_ERR, quantities_rel_L2))

            if rel_tol_fulfilled \
            or  REFERENCE_REFINEMENT >= MAX_IT:
                break

            REFERENCE_REFINEMENT += 1
    else:
        time_mesh_ref    = time_and_log_subprocess(command='python3',  args=['mesh.py','--reference', '--refine', str(REFERENCE_REFINEMENT)], stdout_filename=logfile, stderr_filename=logfile)
        solve_args = ['solve.py','--reference']
        if STEADY_STATE:
            solve_args.append('--steady_state')
        time_solve_ref   = time_and_log_subprocess(command='python3', args=solve_args, stdout_filename=logfile, stderr_filename=logfile)

# mesh and solve mortar+TSA
output_dir_base = 'out/'
files = []
iteration = 0

while(True):
    time_mesh    = time_and_log_subprocess(command='python3',  args=['mesh.py', '--refine', str(iteration)], stdout_filename=logfile, stderr_filename=logfile)
    solve_args = ['solve.py']
    if STEADY_STATE:
        solve_args.append('--steady_state')
    time_solve   = time_and_log_subprocess(command='python3', args=solve_args, stdout_filename=logfile, stderr_filename=logfile)
    ref_name_ext = "_ana.pos" if USE_ANALYTICAL else "_ref.pos"

    # L2
    quantities_ref_L2 = [compute_err_L2(q+ref_name_ext, "", logfile)[-1] for q in QUANTITIES]
    quantities_abs_L2 = [compute_err_L2(q+".pos", q+ref_name_ext, logfile)[-1] for q in QUANTITIES]
    quantities_rel_L2 = [quantities_abs_L2[i] / quantities_ref_L2[i] for i in range(len(QUANTITIES))]

    # L1
    quantities_ref_L1 = [compute_err_L1(q+ref_name_ext, "", logfile)[-1] for q in QUANTITIES]
    quantities_abs_L1 = [compute_err_L1(q+".pos", q+ref_name_ext, logfile)[-1] for q in QUANTITIES]
    quantities_rel_L1 = [quantities_abs_L1[i] / quantities_ref_L1[i] for i in range(len(QUANTITIES))]

    # LInf
    quantities_ref_LInf = [compute_err_LInf(q+ref_name_ext, "", logfile)[-1] for q in QUANTITIES]
    quantities_abs_LInf = [compute_err_LInf(q+".pos", q+ref_name_ext, logfile)[-1] for q in QUANTITIES]
    quantities_rel_LInf = [quantities_abs_LInf[i] / quantities_ref_LInf[i] for i in range(len(QUANTITIES))]

    # Tmax
    quantities_Tmax = [compute_err_LInf(q+".pos", "", logfile)[-1] for q in QUANTITIES]

    # format: display name, reference L2, absolute error, relative error
    labels = ("Reference L2", "Absolute L2", "Relative L2", "Reference L1", "Absolute L1", "Relative L1", "Absolute LInf", "Tmax")

    iteration_data = [[QUANTITIES[i],
        quantities_ref_L2[i], quantities_abs_L2[i], quantities_rel_L2[i],
        quantities_ref_L1[i], quantities_abs_L1[i], quantities_rel_L1[i],
        quantities_abs_LInf[i],
        quantities_Tmax[i]] for i in range(len(QUANTITIES))]
    report = iteration_report(iteration, labels, iteration_data)
    print_and_log(report, 0, logfile)

    # save results
    output_dir = output_dir_base + 'iteration_'+str(iteration)
    files = []
    files += glob('*.pos')
    files += glob('*.msh')
    files += glob('*.log')
    time_and_log_subprocess(command='mkdir', args=['-p', output_dir], stdout_filename=logfile, stderr_filename=logfile)
    time_and_log_subprocess(command='cp', args=['-a'] + files + [output_dir], stdout_filename=logfile, stderr_filename=logfile)

    iteration += 1

    rel_tol_fulfilled = list(filter(lambda x: x < TOL_REL_ERR, quantities_rel_L2))
    abs_tol_fulfilled = list(filter(lambda x: x < TOL_ABS_ERR, quantities_abs_L2))

    if rel_tol_fulfilled \
    or abs_tol_fulfilled \
    or  iteration >= MAX_IT:
        break

# save results
archive_name = 'test_' + str(isodate()) + '.tar.gz'
time_and_log_subprocess(command='tar', args=['czf', archive_name, output_dir_base], stdout_filename=logfile, stderr_filename=logfile)
time_and_log_subprocess(command='rm', args=['-rf', output_dir_base], stdout_filename=logfile, stderr_filename=logfile)

# go back to main dir
os.chdir(main_dir)