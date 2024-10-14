import subprocess
from decimal import Decimal, getcontext
import signal
import sys
#this scrip executes mc by checking statistics
def format_using_decimal(value, precision=4):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)



if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

T=float(sys.argv[1])


N=int(sys.argv[2])


TStr=format_using_decimal(T)
print("TStr="+TStr)
#############################################
#launch mc, i.e., giving initial conditions

launchResult=subprocess.run(["python3", "launch_one_run.py", "./dataAll/dataAllUnitCell"+str(N)+"/T"+TStr+"/run_T"+TStr+".mc.conf"])
print(launchResult.stdout)
if launchResult.returncode!=0:
    print("error in launch one run: "+str(launchResult.returncode))
#############################################


#############################################
#cmake ., make run_mc
targetName="run_mc"
compileErrCode=10
cmake_process=subprocess.Popen(["cmake","."], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = cmake_process.stdout.readline()
    if output == '' and cmake_process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = cmake_process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())

make_process=subprocess.Popen(["make",targetName], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = make_process.stdout.readline()
    if output == '' and make_process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = make_process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())
#############################################



#############################################
#run executable
# Function to handle the termination of the subprocess on SIGINT
def terminate_process(signal_number, frame):
    print("Terminating subprocess...")
    cpp_process.terminate()
    cpp_process.wait()  # Wait for the subprocess to properly terminate
    sys.exit(0)  # Exit the script
cppExecutable="./run_mc"
cpp_process = subprocess.Popen(
    [cppExecutable, "./dataAll/dataAllUnitCell"+str(N)+"/T"+TStr+"/cppIn.txt"],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

# Register the signal handler
signal.signal(signal.SIGINT, terminate_process)

# Read output line by line in real-time
try:
    for line in iter(cpp_process.stdout.readline, ''):
        print(line.strip())

    # Wait for the process to finish and get the final output
    stdout, stderr = cpp_process.communicate()

    # Print any remaining output
    if stdout:
        print(stdout.strip())
    if stderr:
        print(stderr.strip())
except Exception as e:
    print(f"An error occurred: {e}")
finally:
    # Ensure the subprocess is terminated if the script exits unexpectedly
    cpp_process.terminate()
    cpp_process.wait()

#############################################


print("T="+TStr)