import numpy as np
from pandas import DataFrame

from pyzx.circuit import Circuit
from pyzx.circuit.gates import CNOT, gate_types
from pyzx.linalg import Mat2

from pytket import Transform, OpType
from pytket.pyzx import pyzx_to_tk

def get_metrics(circuit):
    metrics = {}
    if isinstance(circuit, Circuit):
        tk_circuit = pyzx_to_tk(circuit)
    else:
        tk_circuit = circuit
    # Make sure everything is in the same gate set!
    Transform.RebaseToQiskit().apply(tk_circuit) # Make all gates into U1, U2, U3, CX
    Transform.RebaseToRzRx().apply(tk_circuit) # Make single qubit gates into Rz, Rx
    metrics["CX depth"] = tk_circuit.depth_by_type(OpType.CX)
    metrics["# CX"] = tk_circuit.n_gates_of_type(OpType.CX)
    metrics["depth"] = tk_circuit.depth()
    metrics["Rz depth"] = tk_circuit.depth_by_type(OpType.Rz)
    metrics["# Rz"] = tk_circuit.n_gates_of_type(OpType.Rz)
    metrics["Rx depth"] = tk_circuit.depth_by_type(OpType.Rx)
    metrics["# Rx"] = tk_circuit.n_gates_of_type(OpType.Rx)
    return DataFrame([metrics])

class CNOT_tracker(Circuit):
    
    def __init__(self, n_qubits, **kwargs):
        super().__init__(n_qubits, **kwargs)
        self.matrix = Mat2(np.identity(n_qubits, dtype=np.int32).tolist())
        self.row_perm = np.arange(n_qubits)
        self.col_perm = np.arange(n_qubits)
        self.n_qubits = n_qubits

    def row_add(self, q0, q1):
        self.add_gate("CNOT", q0, q1)
        self.matrix.row_add(q0, q1)

    def add_gate(self, gate, *args, **kwargs):
        if isinstance(gate, CNOT):
            self.row_add(gate.control, gate.target)
        else:
            super().add_gate(gate, *args, **kwargs)

    def col_add(self, q0, q1):
        self.prepend_gate("CNOT", q1, q0)
        self.matrix.col_add(q0, q1)

    def prepend_gate(self, gate, *args, **kwargs):
        """Adds a gate to the circuit. ``gate`` can either be 
        an instance of a :class:`Gate`, or it can be the name of a gate,
        in which case additional arguments should be given.

        Example::

            circuit.add_gate("CNOT", 1, 4) # adds a CNOT gate with control 1 and target 4
            circuit.add_gate("ZPhase", 2, phase=Fraction(3,4)) # Adds a ZPhase gate on qubit 2 with phase 3/4
        """
        if isinstance(gate, str):
            gate_class = gate_types[gate]
            gate = gate_class(*args, **kwargs)
        self.gates.insert(0, gate)

    def to_qasm(self):
        qasm = super().to_qasm()
        initial_perm = "// Initial wiring: " + str(self.row_perm.tolist())
        end_perm = "// Resulting wiring: " + str(self.col_perm.tolist())
        return '\n'.join([initial_perm, end_perm, qasm])

    @staticmethod
    def from_circuit(circuit):
        new_circuit = CNOT_tracker(circuit.qubits, name=circuit.name)
        new_circuit.gates = circuit.gates
        new_circuit.update_matrix()
        return new_circuit

    def update_matrix(self):
        self.matrix = Mat2(np.identity(self.n_qubits))
        for gate in self.gates:
            if hasattr(gate, "name") and gate.name == "CNOT":
                self.matrix.row_add(gate.control, gate.target)
            else:
                print("Warning: CNOT tracker can only be used for circuits with only CNOT gates!")

    @staticmethod
    def from_qasm_file(fname):
        circuit = Circuit.from_qasm_file(fname)
        return CNOT_tracker.from_circuit(circuit)