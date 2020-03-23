from pytket import route, Architecture, Transform, OpType

from pyzx.circuit import Circuit
from pyzx.circuit.gates import ZPhase, CNOT

def route_tket(circuit, architecture):
    outcirc = route(circuit, architecture)
    Transform.OptimisePostRouting().apply(outcirc)
    return outcirc

def get_tk_architecture(architecture):
    coupling_graph = [e for e in architecture.graph.edges()]
    return Architecture(coupling_graph)
    