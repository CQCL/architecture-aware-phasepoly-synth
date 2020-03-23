from pandas import DataFrame, concat
import numpy as np 
from fractions import Fraction

from pytket import PauliExpBox, Pauli, Transform
from pytket.pyzx import tk_to_pyzx
from pytket import Circuit as TketCircuit

from pyzx.circuit import Circuit
from pyzx.circuit.gates import ZPhase, HAD, CNOT
from pyzx.linalg import Mat2

from steiner import steiner_reduce_column
from steiner import rec_steiner_gauss as gauss
from architecture import create_architecture, FULLY_CONNECTED
from utils import CNOT_tracker
        
def route_phase_poly(circuit, architecture, do_matroid="gray", **kwargs):
    if not isinstance(circuit, PhasePoly):
        phase_poly = PhasePoly.fromCircuit(circuit)
    else:
        phase_poly = circuit
    new_circuit = None
    if do_matroid == "arianne":
        new_circuit = phase_poly.Ariannes_synth(architecture, **kwargs)[0]
    elif do_matroid == "gray":
        new_circuit = phase_poly.gray_synth(architecture, **kwargs)[0]
    return new_circuit

def make_random_phase_poly(n_qubits, n_gadgets, return_circuit=False):
    parities = set()
    if n_gadgets > 2**n_qubits:
        n_gadgets=n_qubits^3
    if n_qubits < 26:
        for integer in np.random.choice(2**n_qubits-1, replace=False, size=n_gadgets):
            parities.add(integer+1)
    elif n_qubits < 64:
        while len(parities) < n_gadgets:
            parities.add(np.random.randint(1, 2**n_qubits))
    else:
        while len(parities) < n_gadgets:
            parities.add("".join(np.random.choice(["0", "1"], n_qubits, replace=True)))
    if n_qubits < 64:
        parities = [("{0:{fill}"+str(n_qubits)+"b}").format(integer, fill='0', align='right') for integer in parities]
    zphase_dict = {"".join([str(int(i)) for i in p]):Fraction(1,4) for p in parities}
    out_parities = mat22partition(Mat2.id(n_qubits))
    phase_poly = PhasePoly(zphase_dict, out_parities)
    if return_circuit:
        return route_phase_poly(phase_poly, create_architecture(FULLY_CONNECTED, n_qubits=n_qubits), do_matroid="arianne")
    return phase_poly

class PhasePoly():

    def __init__(self, zphase_dict, out_parities):
        self.zphases = zphase_dict
        self.out_par = out_parities
        self.n_qubits = len(out_parities[0])
        self.all_parities = list(zphase_dict.keys())

    @staticmethod
    def fromCircuit(circuit, initial_qubit_placement=None, final_qubit_placement=None):
        if isinstance(circuit, TketCircuit):
            circuit = tk_to_pyzx(circuit)
        zphases = {}
        current_parities = mat22partition(Mat2.id(circuit.qubits))
        if initial_qubit_placement is not None:
            current_parities = ["".join([row[i] for i in initial_qubit_placement]) for row in current_parities]
        for gate in circuit.gates:
            parity = current_parities[gate.target]
            if gate.name in ["CNOT", "CX"]:
                # Update current_parities
                control = current_parities[gate.control]
                current_parities[gate.target] = "".join([str((int(i)+int(j))%2) for i,j in zip(control, parity)])
            elif isinstance(gate, ZPhase):
                # Add the T rotation to the phases
                if parity in zphases:
                    zphases[parity] += gate.phase
                else: 
                    zphases[parity] = gate.phase
            else:
                print("Gate not supported!", gate.name)
        def clamp(phase):
            new_phase = phase%2
            if new_phase > 1:
                return new_phase -2
            return new_phase
        zphases = {par:clamp(r) for par, r in zphases.items() if clamp(r) != 0}
        if final_qubit_placement is not None:
            current_parities = [ current_parities[i] for i in final_qubit_placement]
        return PhasePoly(zphases, current_parities)

    def to_tket(self):
        if self.out_par != mat22partition(Mat2.id(self.n_qubits)):
            print(self.out_par)
            raise NotImplementedError("The linear transformation part of the phase polynomial cannot yet be transformed into a tket circuit")
        circuit = TketCircuit(self.n_qubits)
        for parity, phase in self.zphases.items():
            qubits = [i for i,s in enumerate(parity) if s == '1']
            circuit.add_pauliexpbox(PauliExpBox([Pauli.Z]*len(qubits), phase), qubits)
        # TODO add the linear combination part with CNOTs
        Transform.DecomposeBoxes().apply(circuit)
        return circuit

    def gray_synth(self, architecture, full=True, **kwargs):
        if architecture is None:
            architecture = create_architecture(FULLY_CONNECTED, n_qubits=len(self.out_par[0]))
        n_qubits = architecture.n_qubits
        # Obtain the parities
        parities_to_reach = self.all_parities
        # Make a matrix from the parities
        matrix = Mat2([[int(parity[i]) for parity in parities_to_reach] for i in range(n_qubits)] )
        circuit = CNOT_tracker(n_qubits)
        # Make a stack - aka use the python stack >^.^<
        def recurse(cols_to_use, qubits_to_use, phase_qubit): # Arguments from the original paper
            # Check for finished columns
            cols_to_use = self._check_columns(matrix, circuit, cols_to_use, parities_to_reach)
            if cols_to_use != [] and qubits_to_use != []:   
                # Find all qubits (rows) with only 1s on the allowed parities (cols_to_use) 
                qubits = [i for i in qubits_to_use if sum([matrix.data[i][j] for j in cols_to_use]) == len(cols_to_use)]
                if len(qubits) > 1 and phase_qubit is not None:
                    # Pick the column with the most 1s to extrac the steiner tree with
                    column = max(cols_to_use, key=lambda c: sum([row[c] for row in matrix.data]))
                    col = [1 if i in qubits else 0 for i in range(n_qubits)]
                    cnots =  list(steiner_reduce_column(architecture, col, phase_qubit, set(qubits + [phase_qubit]), [i for i in range(n_qubits)], [], upper=True))
                    # For each returned CNOT:
                    for target, control in cnots:
                        # Place the CNOT on the circuit
                        circuit.add_gate(CNOT(control, target))
                        # Adjust the matrix accordingly - reversed elementary row operations
                        matrix.row_add(target, control)
                        # Keep track of the parities in the circuit - normal elementary row operations
                        cols_to_use = self._check_columns(matrix, circuit, cols_to_use, parities_to_reach)   
                # After placing the cnots do recursion
                if len(cols_to_use) > 0:
                    # Choose a row to split on
                    if len(cols_to_use) == 1: 
                        # Hack because it runs into infinite recursion otherwise, because we do not remove the chosen_row from qubits_to_use
                        # We do not remove the qubit, because it will skip columns 
                        # It skips columns because it runs out of qubits before all columns are finished for some reason - could not find the bug
                        qubits = [i for i in qubits_to_use if sum([matrix.data[i][j] for j in cols_to_use]) == 1]
                    else:
                        # Ignore rows that are currently all 0s or all 1s
                        qubits = [i for i in qubits_to_use if sum([matrix.data[i][j] for j in cols_to_use]) not in [0, len(cols_to_use)]]
                    # Pick the qubit where the recursion split will be most skewed.
                    chosen_row = max(qubits, key=lambda q: max([len([col for col in cols_to_use if matrix.data[q][col] == i]) for i in [1,0]], default=-1))
                    # Split the column into 1s and 0s in that row
                    cols1 = [col for col in cols_to_use if matrix.data[chosen_row][col] == 1]
                    cols0 = [col for col in cols_to_use if matrix.data[chosen_row][col] == 0]
                    #qubits_to_use = [i for i in qubits_to_use if i != chosen_row]
                    recurse(cols0, qubits_to_use, phase_qubit)
                    recurse(cols1, qubits_to_use, phase_qubit if phase_qubit is not None else chosen_row)
        # Put the base case into the python stack
        recurse([i for i in range(len(parities_to_reach))], [i for i in range(n_qubits)], None)
        if full:
            # Calculate the final parity that needs to be added from the circuit and self.out_par
           self._obtain_final_parities(circuit, architecture, **kwargs)
        # Return the circuit
        circuit.n_gadgets = len(self.zphases.keys())
        return circuit, [i for i in range(architecture.n_qubits)], [i for i in range(architecture.n_qubits)]

    def Ariannes_synth(self, architecture, full=True, depth_aware=False, **kwargs):
        kwargs["full_reduce"] = True
        #print("------------ START! -----------------------")
        if architecture is None:
            architecture = create_architecture(FULLY_CONNECTED, n_qubits=len(self.out_par[0]))
        n_qubits = architecture.n_qubits
        # Obtain the parities
        parities_to_reach = self.all_parities
        # Make a matrix from the parities
        matrix = Mat2([[int(parity[i]) for parity in parities_to_reach] for i in range(architecture.n_qubits)] )
        circuit = CNOT_tracker(architecture.n_qubits)
        cols_to_reach = self._check_columns(matrix, circuit, [i for i in range(len(parities_to_reach))], parities_to_reach)
        self.prev_rows = None
        def place_cnot(control, target):
            # Place the CNOT on the circuit
            circuit.add_gate(CNOT(control, target))
            # Adjust the matrix accordingly - reversed elementary row operations
            matrix.row_add(target, control) 
        def base_recurse(cols_to_use, qubits_to_use):
            if qubits_to_use != [] and cols_to_use != []:
                # Select edge qubits
                qubits = architecture.non_cutting_vertices(qubits_to_use) 
                # Pick the qubit where the recursion split will be most skewed.
                chosen_row = max(qubits, key=lambda q: max([len([col for col in cols_to_use if matrix.data[q][col] == i]) for i in [1, 0]], default=-1))
                # Split the column into 1s and 0s in that row
                cols1 = [col for col in cols_to_use if matrix.data[chosen_row][col] == 1]
                cols0 = [col for col in cols_to_use if matrix.data[chosen_row][col] == 0]
                base_recurse(cols0, [q for q in qubits_to_use if q != chosen_row])
                one_recurse(cols1, qubits_to_use, chosen_row)
        def one_recurse(cols_to_use, qubits_to_use, qubit):
            if cols_to_use != []:
                neighbors = [q for q in architecture.get_neighboring_qubits(qubit) if q in qubits_to_use ]
                if neighbors == []:
                    print(qubits_to_use, qubit)
                    print(matrix)
                    exit(42)
                chosen_neighbor = max(neighbors, key=lambda q: len([col for col in cols_to_use if matrix.data[q][col] == 1]))
                # Place CNOTs if you still need to extract columns
                if sum([matrix.data[chosen_neighbor][c] for c in cols_to_use]) != 0: # Check if adding the cnot is useful
                    place_cnot(qubit, chosen_neighbor)
                    # Might have changed the matrix.
                    cols_to_use = self._check_columns(matrix, circuit, cols_to_use, parities_to_reach)  
                else: # Will never change the matrix
                    place_cnot(chosen_neighbor, qubit)
                    place_cnot(qubit, chosen_neighbor)
                    # Since the neighbor was all zeros, this is effectively a swap and no columns need to be checked.
                # Split the column into 1s and 0s in that row
                cols0 = [col for col in cols_to_use if matrix.data[qubit][col] == 0]
                cols1 = [col for col in cols_to_use if matrix.data[qubit][col] == 1]
                base_recurse(cols0, [q for q in qubits_to_use if q != qubit])
                one_recurse(cols1, qubits_to_use, qubit)
        base_recurse(cols_to_reach, [i for i in range(n_qubits)])
        if full:
            # Calculate the final parity that needs to be added from the circuit and self.out_par
           self._obtain_final_parities(circuit, architecture, **kwargs)
        # Return the circuit
        circuit.n_gadgets = len(self.zphases.keys())
        return circuit, [i for i in range(architecture.n_qubits)], [i for i in range(architecture.n_qubits)]

    def _check_columns(self, matrix, circuit, columns, parities_to_reach):
        # Check for finished columns
        for col in [c for c in columns]:
            column_data = [row[col] for row in matrix.data]
            if sum(column_data) == 1:
                # Add phase gates where needed
                qubit = column_data.index(1)
                circuit.add_gate(ZPhase(qubit, self.zphases[parities_to_reach[col]]))
                # Remove columns from the matrix if the corresponding parity was obtained
                columns.remove(col)
        return columns 

    def _obtain_final_parities(self, circuit, architecture, **kwargs):
        # Calculate the final parity that needs to be added from the circuit and self.out_par
        current_parities = circuit.matrix
        output_parities = Mat2([[int(v) for v in row] for row in self.out_par])
        last_parities = output_parities*current_parities.inverse()
        # Do steiner-gauss to calculate necessary CNOTs and add those to the circuit.
        cnots = CNOT_tracker(architecture.n_qubits)
        gauss(last_parities, architecture, y=cnots, **kwargs) 
        for cnot in cnots.gates:
            circuit.add_gate(cnot)
            
def partition2mat2(partition):
    return Mat2([[int(i) for i in parity] for parity in partition])

def mat22partition( m):
    return ["".join(str(i) for i in parity) for parity in m.data]
