import subprocess, tempfile, datetime

from pytket import circuit_from_qasm, Transform
from pytket.pyzx import pyzx_to_tk

from tket_router import route_tket, get_tk_architecture
from phase_poly import route_phase_poly, PhasePoly

PROPOSED = "proposed"
NASH = "nash"
TKET = "tket"
STAQ = "staq"

staq_device_map = {
    "rigetti_16q_aspen": "aspen-4",
    "9q-square":"square",
    "ibm_q20_tokyo":"tokyo",
    "9q-fully_connected":"fullcon",
    "rigetti_8q_agave":"agave",
    "ibmq_singapore":"singapore"
}

def route_nash(circuit, architecture):
    return route_phase_poly(circuit, architecture, do_matroid="gray")

def route_proposed(circuit, architecture):
    return route_phase_poly(circuit, architecture, do_matroid="arianne")

def route_with_tket(circuit, tket_architecture):
    circuit = PhasePoly.fromCircuit(circuit).to_tket()
    return route_tket(circuit,tket_architecture)

def route_staq(file, staq_arch_name, mapped=False, return_circuit=True):
    if mapped:
        base_command = "./staq -M steiner -m -d"
    else:
        base_command = "./staq -M steiner --disable-layout-optimization -m -d"
    s = subprocess.check_output(" ".join([base_command, staq_arch_name, file]), shell=True)
    if return_circuit:
        s = s.decode("utf-8")
        s = s.replace("U", "u3")
        s = s.replace("CX", "cx")
        fp = tempfile.NamedTemporaryFile(suffix=".qasm") # Keep the qasm file in memory, but pretend it exists on disk
        fp.write(s.encode('utf-8'))
        fp.seek(0) #Place the pointer back at the start of the file
        try:
            c = circuit_from_qasm(fp.name)
        except TypeError as e:
            print(s)
            raise e
        fp.close()
        return c
    else:
        return s

def route_with(filename, architecture, mode, mapped=False):
    circuit = circuit_from_qasm(filename)
    if mode == PROPOSED:
        start = datetime.datetime.now()
        c = route_proposed(circuit, architecture)
        time = datetime.datetime.now() - start
    elif mode == NASH:
        start = datetime.datetime.now()
        c = route_nash(circuit, architecture)
        time = datetime.datetime.now() - start
    elif mode == TKET:
        tk_arch = get_tk_architecture(architecture)
        start = datetime.datetime.now()
        c = route_with_tket(circuit, tk_arch)
        time = datetime.datetime.now() - start
    elif mode == STAQ:
        staq_device = staq_device_map[architecture.name]
        start = datetime.datetime.now()
        s = route_staq(filename, staq_device, mapped=mapped, return_circuit=False)
        time = datetime.datetime.now() - start
        
        s = s.decode("utf-8")
        s = s.replace("U", "u3")
        s = s.replace("CX", "cx")
        fp = tempfile.NamedTemporaryFile(suffix=".qasm") # Keep the qasm file in memory, but pretend it exists on disk
        fp.write(s.encode('utf-8'))
        fp.seek(0) #Place the pointer back at the start of the file
        try:
            c = circuit_from_qasm(fp.name)
        except TypeError as e:
            print(s)
            raise e
        fp.close()
    return c, time