import sys
from pyzx.graph.graph import  Graph
import numpy as np
import networkx as nx

SQUARE = "square"
LINE = "line"
FULLY_CONNECTED = "fully_connected"
CIRCLE = "circle"
IBM_QX2 = "ibm_qx2"
IBM_QX3 = "ibm_qx3"
IBM_QX4 = "ibm_qx4"
IBM_QX5 = "ibm_qx5"
IBM_Q20_TOKYO = "ibm_q20_tokyo"
RIGETTI_19Q_ACORN = "rigetti_19q_acorn"
RIGETTI_16Q_ASPEN = "rigetti_16q_aspen"
RIGETTI_8Q_AGAVE = "rigetti_8q_agave"
REC_ARCH = "recursive_architecture"
SYCAMORE_LIKE = "sycamore_like"
IBMQ_POUGHKEEPSIE = "ibmq_poughkeepsie"
IBMQ_SINGAPORE = "ibmq_singapore"
IBMQ_BOEBLINGEN  = "ibmq_boeblingen"
GOOGLE_SYCAMORE = "google_sycamore"
IBM_ROCHESTER = "ibm_rochester"
DENSITY = "dynamic_density"

architectures = [SQUARE, CIRCLE, FULLY_CONNECTED, LINE, DENSITY, IBM_QX4, IBM_QX2, IBM_QX3, 
                IBM_QX5, IBM_Q20_TOKYO, RIGETTI_8Q_AGAVE, RIGETTI_16Q_ASPEN, RIGETTI_19Q_ACORN, 
                REC_ARCH, SYCAMORE_LIKE, IBMQ_POUGHKEEPSIE, IBMQ_BOEBLINGEN, IBMQ_SINGAPORE, 
                GOOGLE_SYCAMORE, IBM_ROCHESTER]
dynamic_size_architectures = [FULLY_CONNECTED, LINE, CIRCLE, SQUARE, DENSITY]

debug = False

class Architecture():
    def __init__(self, name, coupling_graph=None, coupling_matrix=None, backend=None, qubit_map = None, reduce_order=None, **kwargs):
        """
        Class that represents the architecture of the qubits to be taken into account when routing.

        :param coupling_graph: a PyZX Graph representing the architecture, optional 
        :param coupling_matrix: a 2D numpy array representing the adjacency of the qubits, from which the Graph is created, optional
        :param backend: The PyZX Graph backend to be used when creating it from the adjacency matrix, optional
        """
        self.name = name
        if coupling_graph is None:
            self.graph = Graph(backend=backend)
            self._double_edged = True
        else:
            self.graph = coupling_graph
            self._double_edged = False

        if coupling_matrix is not None:
            # build the architecture graph
            n = coupling_matrix.shape[0]
            self.vertices = self.graph.add_vertices(n)
            edges = [(self.vertices[row], self.vertices[col]) for row in range(n) for col in range(n) if
                     coupling_matrix[row, col] == 1]
            self.graph.add_edges(edges)
        else:
            self.vertices = [v for v in self.graph.vertices()]
        self.distances = None
        if qubit_map is not None:
            self.qubit_map = qubit_map
        else:
            self.qubit_map = [i for i, v in enumerate(self.vertices)]
        self.n_qubits = len(self.vertices)
        self._non_cutting_vertices = None
        #self.pre_calc_non_cutting_vertices()
        self.reduce_order = reduce_order if reduce_order is not None else [i-1 for i in range(self.n_qubits, 0,-1)]
        self.density = self._density(self.graph.nedges, self.n_qubits)

    def pre_calc_distances(self):
        self.distances = {"upper": [self.floyd_warshall(until, upper=True) for until, v in enumerate(self.vertices)],
                          "full": [self.floyd_warshall(until, upper=False) for until, v in enumerate(self.vertices)]}

    def _density(self, edges, qubits):
        density = edges/(qubits*(qubits-1))
        if self._double_edged:
            return density
        return 2*density

    def pre_calc_non_cutting_vertices(self):
        # TODO implement this and uncomment line in constructor.
        #raise NotImplementedError("pre calculation non cutting vertices")
        qubits = [i for i in range(self.n_qubits)]
        def collect_non_cutting(qubits):
            if qubits == []:
                return []
            vertices = [self.vertices[q] for q in qubits]
            is_cutting = self._is_cutting(vertices)
            non_cutting = [q for i,q  in enumerate(qubits) if not is_cutting[i] ]
            all_non_cutting = [(qubits, non_cutting)]
            for qubit in non_cutting:
                subgraph = [q for q in qubits if q != qubit]
                all_non_cutting += collect_non_cutting(subgraph)
            return all_non_cutting
        self._non_cutting_vertices = collect_non_cutting(qubits)
    
    def non_cutting_vertices(self, subgraph, pre_calc=False):
        # Find the precalculated non-cutting vertices for this subgraph
        if pre_calc:
            if self._non_cutting_vertices is None:
                self.pre_calc_non_cutting_vertices()
            if subgraph == []:
                return []
            subgraphs, non_cutting = zip(*self._non_cutting_vertices)
            index = subgraphs.index(subgraph)
            return non_cutting[index]
        else:
            if self._non_cutting_vertices is None:
                self._non_cutting_vertices = {}
            cur_dict = self._non_cutting_vertices
            for q in sorted([v for i, v in enumerate(self.vertices) if i not in subgraph]):
                if q not in cur_dict.keys():
                    cur_dict[q] = {}
                cur_dict = cur_dict[q]
            if "non_cutting" not in cur_dict.keys():
                cur_dict["non_cutting"] = [subgraph[i] for i, cutting in enumerate(self._is_cutting(qubits=[self.vertices[i] for i in subgraph])) if not cutting]
            #print(self._non_cutting_vertices)
            return cur_dict["non_cutting"]

    def _is_cutting(self, qubits=None):
        # algortihm from https://courses.cs.washington.edu/courses/cse421/04su/slides/artic.pdf and https://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/ 
        # Code written myself.
        if qubits is None:
            qubits = self.vertices
        n_qubits = len(qubits)
        discovery_times = [-1]*n_qubits
        lows = [None]*n_qubits
        index_lookup = {self.vertices[v]:i for i, v in enumerate(qubits)}
        self.dfs_counter = 0
        edges = [e for e in self.graph.edges() if e[0] in qubits and e[1] in qubits]
        edges += [(v2, v1) for v1, v2 in edges]
        cutting = [False] * n_qubits
        parent = [-1] * n_qubits
        def dfs(vertex,):
            v = index_lookup[vertex]
            self.dfs_counter += 1
            discovery_times[v] = self.dfs_counter
            lows[v] = discovery_times[v]
            children = 0
            for edge in [e for e in edges if e[0] == vertex]:
                vertex2 = edge[1]
                v2 = index_lookup[vertex2]
                if discovery_times[v2] == -1: # Not visited yet
                    parent[v2] = v
                    children += 1
                    dfs(vertex2)
                    lows[v] = min(lows[v], lows[v2])
                    if parent[v] == -1 and children > 1:
                        cutting[v] = True
                    if parent[v] != -1 and lows[v2] >= discovery_times[v]:
                        cutting[v] = True
                elif v2 != parent[v]:
                    lows[v] = min(lows[v], discovery_times[v2])
        for vertex in qubits:
            v = index_lookup[vertex]
            if discovery_times[v] == -1:
                dfs(vertex)
        del self.dfs_counter
        return cutting

    def get_neighboring_qubits(self, qubit):
        neighboring_qubits = set([edge[1] for edge in self.graph.edges() if edge[0] == self.vertices[qubit] ] + [edge[0] for edge in self.graph.edges() if edge[1] == self.vertices[qubit]])
        return [self.vertices.index(q) for q in neighboring_qubits]

    def visualize(self, filename=None):
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')
        g = nx.Graph()
        g.add_nodes_from(self.vertices)
        g.add_edges_from(self.graph.edges())
        nx.draw(g, with_labels=True, font_weight='bold')
        if filename is None:
            filename = self.name + ".png"
        plt.savefig(filename)

    def floyd_warshall(self, exclude_excl, upper=True):
        """
        Implementation of the Floyd-Warshall algorithm to calculate the all-pair distances in a given graph

        :param exclude_excl: index up to which qubit should be excluded from the distances
        :param upper: whether use bidirectional edges or only ordered edges (src, tgt) such that src > tgt, default True
        :return: a dict with for each pair of qubits in the graph, a tuple with their distance and the corresponding shortest path
        """
        # https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
        distances = {}
        vertices = self.vertices[exclude_excl:] if upper else self.vertices[:exclude_excl + 1]
        for edge in self.graph.edges():
            src, tgt = self.graph.edge_st(edge)
            if src in vertices and tgt in vertices:
                if upper:
                    distances[(src, tgt)] = (1, [(src, tgt)])
                    distances[(tgt, src)] = (1, [(tgt, src)])
                elif src > tgt:
                    distances[(src, tgt)] = (1, [(src, tgt)])
                else:
                    distances[(tgt, src)] = (1, [(tgt, src)])
        for v in vertices:
            distances[(v, v)] = (0, [])
        for i, v0 in enumerate(vertices):
            for j, v1 in enumerate(vertices if upper else vertices[:i + 1]):
                for v2 in vertices if upper else vertices[: i + j + 1]:
                    if (v0, v1) in distances.keys():
                        if (v1, v2) in distances.keys():
                            if (v0, v2) not in distances.keys() or distances[(v0, v2)][0] > distances[(v0, v1)][0] + \
                                    distances[(v1, v2)][0]:
                                distances[(v0, v2)] = (distances[(v0, v1)][0] + distances[(v1, v2)][0],
                                                       distances[(v0, v1)][1] + distances[(v1, v2)][1])
                                if upper:
                                    distances[(v2, v0)] = (distances[(v0, v1)][0] + distances[(v1, v2)][0],
                                                       distances[(v2, v1)][1] + distances[(v1, v0)][1])
        return distances

    def shortest_path(self, start, end, subgraph=None):
        #print("shortest path")
        if subgraph is None:
            nodes = self.vertices
        else: 
            nodes = [self.vertices[n] for n in subgraph]
        start = self.vertices[start]
        end = self.vertices[end]

        queue = [(start, [start])] 
        visited = [start]
        debug_trace = [queue[0]]
  
        while queue != []: 
            
            node, path = queue.pop() 
            if node == end:
                #print(start, end, subgraph)
                #self.visualize("debug.png")
                #print(path)
                return path
            edges = [edge for edge in self.graph.edges() if node in edge]
            neighbors = [n for edge in edges for n in edge if n != node and n in nodes]
            #print(node, neighbors)
            for node in neighbors: 
                if node not in visited: 
                    queue.append((node, path + [node]))
                    visited.append(node)
                    debug_trace.append(queue[-1])
        print(start, end, subgraph)
        print(*self.graph.edges())
        print("Trace:", *debug_trace, sep="\n")
        self.visualize("debug.png")
        print(self.name)
        return None

    def steiner_tree(self, start, nodes, upper=True):
        """
        Approximates the steiner tree given the architecture, a root qubit and the other qubits that should be present.
        This is done using the pre-calculated all-pairs shortest distance and Prim's algorithm for creating a minimum spanning tree
        :param start: The index of the root qubit to be used
        :param nodes: The indices of the other qubits that should be present in the steiner tree
        :param upper: Whether the steiner tree is used for creating an upper triangular matrix or a full reduction.
        :yields: First yields all edges from the tree top-to-bottom, finished with None, then yields all edges from the tree bottom-up, finished with None.
        """
        # Approximated by calculating the all-pairs shortest paths and then solving the mininum spanning tree over the subset of vertices and their respective shortest paths.
        # https://en.wikipedia.org/wiki/Steiner_tree_problem#Approximating_the_Steiner_tree

        # The all-pairs shortest paths are pre-calculated and the mimimum spanning tree is solved with Prim's algorithm
        # https://en.wikipedia.org/wiki/Prim%27s_algorithm

        # returns an iterator that walks the steiner tree, yielding (adj_node, leaf) pairs. If the walk is finished, it yields None
        if self.distances is None:
            self.pre_calc_distances()
        state = [start, [n for n in nodes]]
        root = start
        # TODO deal with qubit mapping
        vertices = [root]
        edges = []
        debug and print(root, upper, nodes)
        distances = self.distances["upper"][root] if upper else self.distances["full"][root]
        steiner_pnts = []
        while nodes != []:
            options = [(node, v, *distances[(v, node)]) for node in nodes for v in (vertices + steiner_pnts) if
                       (v, node) in distances.keys()]
            best_option = min(options, key=lambda x: x[2])
            debug and print("Adding to tree: vertex ", best_option[0], "Edges ", best_option[3])
            vertices.append(best_option[0])
            edges.extend(best_option[3])
            steiner = [v for edge in best_option[3] for v in edge if v not in vertices]
            debug and print(steiner)
            steiner_pnts.extend(steiner)
            nodes.remove(best_option[0])
        edges = set(edges)  # remove duplicates
        if debug:
            print("edges:", edges)
            print("nodes:", vertices)
            print("steiner points:", steiner_pnts)
        # First go through the tree to find and remove zeros
        state += [[e for e in edges], [v for v in vertices], [s for s in steiner_pnts]]
        vs = {root}
        n_edges = len(edges)
        yielded_edges = set()
        debug_count = 0
        yield_count = 0
        warning = 0
        while len(yielded_edges) < n_edges:
            es = [e for e in edges for v in vs if e[0] == v]
            old_vs = [v for v in vs]
            yielded = False
            for edge in es:
                yield edge
                vs.add(edge[1])
                if edge in yielded_edges:
                    print("DOUBLE yielding! - should not be possible!")
                yielded_edges.add(edge)
                yielded = True
                yield_count += 1
            [vs.remove(v) for v in old_vs]
            if not yielded:
                debug and print("leaf!")
                debug_count += 1
                if debug_count > len(vertices):
                    print("infinite loop!", warning)
                    warning += 1
            if yield_count > len(edges):
                print("Yielded more edges than existing... This should not be possible!", warning)
                warning += 1
            if warning > 5:
                print(state, yielded_edges)
                # input("note it down")
                break
        yield None
        # Walk the tree bottom up to remove all ones.
        yield_count = 0
        while len(edges) > 0:
            # find leaf nodes:
            debug and print(vertices, steiner_pnts, edges)
            vs_to_consider = [vertex for vertex in vertices if vertex not in [e0 for e0, e1 in edges]] + \
                             [vertex for vertex in steiner_pnts if vertex not in [e0 for e0, e1 in edges]]
            yielded = False
            for v in vs_to_consider:
                # Get the edge that is connected to this leaf node
                for edge in [e for e in edges if e[1] == v]:
                    yield edge
                    edges.remove(edge)
                    yielded = True
                    yield_count += 1
                    # yield map(lambda i: self.qubit_map[i], edge)
            if not yielded:
                print("Infinite loop!", warning)
                warning += 1
            if yield_count > n_edges:
                print("Yielded more edges than existing again... This should not be possible!!", warning)
                warning += 1
            if warning > 10:
                print(state, edges, yield_count)
                # input("Note it down!")
                break
        yield None

    def rec_steiner_tree(self, start, nodes, usable_nodes, rec_nodes, upper=True):
        # Builds the steiner tree with start as root, contains at least nodes and at most useable_nodes
        debug_trace = {
            "start":start,
            "nodes": nodes,
            "usable nodes": usable_nodes,
            "rec nodes": rec_nodes,
            "upper": upper
        }
        start = self.qubit_map[start]
        usable_nodes = [self.qubit_map[i] for i in usable_nodes]
        nodes = [self.qubit_map[i] for i in nodes]
        rec_nodes = [self.qubit_map[i] for i in rec_nodes]
        # Calculate all-pairs shortest path
        distances = {}
        vertices = [self.vertices[i] for i in usable_nodes]
        for edge in self.graph.edges():
            src, tgt = self.graph.edge_st(edge)
            if src in vertices and tgt in vertices:
                if upper or (src in rec_nodes and tgt in rec_nodes):
                    distances[(src, tgt)] = (1, [(src, tgt)])
                    distances[(tgt, src)] = (1, [(tgt, src)])
                elif src > tgt:
                    distances[(src, tgt)] = (1, [(src, tgt)])
                else:
                    distances[(tgt, src)] = (1, [(tgt, src)])
        for v in vertices:
            distances[(v, v)] = (0, [])
        for i, v0 in enumerate(vertices+vertices):
            for j, v1 in enumerate(vertices): 
                for v2 in vertices: 
                    if (v0, v1) in distances.keys():
                        if (v1, v2) in distances.keys():
                            if (v0, v2) not in distances.keys() \
                                    or distances[(v0, v2)][0] > distances[(v0, v1)][0] + distances[(v1, v2)][0]:
                                distances[(v0, v2)] = (distances[(v0, v1)][0] + distances[(v1, v2)][0],
                                                        distances[(v0, v1)][1] + distances[(v1, v2)][1])
                                if upper:
                                    distances[(v2, v0)] = (distances[(v0, v1)][0] + distances[(v1, v2)][0],
                                                        distances[(v2, v1)][1] + distances[(v1, v0)][1])
        # Build the spanning tree of shortest paths with root start, containing at least nodes
        vertices = [start]
        edges = []
        steiner_pnts = []
        while nodes != []:
            options = [(node, v, *distances[(v, node)]) for node in nodes for v in (vertices + steiner_pnts) if
                        (v, node) in distances.keys()]
            if options == []:
                print(nodes)
                print(steiner_pnts, vertices)
                print(distances)
                print(usable_nodes)
                print(rec_nodes)
                print(self.reduce_order)
                print(debug_trace)
                self.visualize("debug.png")
                print(self.name)
            best_option = min(options, key=lambda x: x[2])
            vertices.append(best_option[0])
            edges += best_option[3]
            steiner = [v for edge in best_option[3] for v in edge if v not in vertices]
            steiner_pnts += steiner
            nodes.remove(best_option[0])
        edges = list(set(edges)) #removes duplicates

        vs = {start} # Start with the root
        n_edges = len(edges)
        yielded_edges = set()
        while len(yielded_edges) < n_edges:
            es = [e for e in edges for v in vs if e[0] == v] # Find all vertices connected to previously yielded vertices
            old_vs = [v for v in vs]
            for edge in es: # yield the corresponding edges.
                yield (self.qubit_map.index(edge[0]), self.qubit_map.index(edge[1]))
                vs.add(edge[1])
                yielded_edges.add(edge)
            [vs.remove(v) for v in old_vs]
        yield None # Signal next phase
        # Walk the tree bottom up to remove all ones.
        while len(edges) > 0:
            # find leaf nodes:
            vs_to_consider = [vertex for vertex in vertices if vertex not in [e0 for e0, e1 in edges]] + \
                                [vertex for vertex in steiner_pnts if vertex not in [e0 for e0, e1 in edges]]
            for v in vs_to_consider:
                # Get the edge that is connected to this leaf node
                for edge in [e for e in edges if e[1] == v]:
                    yield (self.qubit_map.index(edge[0]), self.qubit_map.index(edge[1])) # yield it
                    edges.remove(edge) # Remove it from the steiner tree
        yield None # Signal done

    def transpose(self):
        # TODO make a transposed copy of self
        qubit_map = list(reversed(self.qubit_map))
        arch = Architecture(self.name + "_transpose", coupling_graph=self.graph, qubit_map=qubit_map)
        return arch

def dynamic_size_architecture_name(base_name, n_qubits):
    return str(n_qubits) + "q-" + base_name

def connect_vertices_in_line(vertices):
    return [(vertices[i], vertices[i+1]) for i in range(len(vertices)-1)]

def connect_vertices_as_grid(width, height, vertices):
    if len(vertices) != width * height:
        raise KeyError("To make a grid, you need vertices exactly equal to width*height, but got %d=%d*%d." % (len(vertices), width, height))
    edges = connect_vertices_in_line(vertices)
    horizontal_lines = [vertices[i*width: (i+1)*width] for i in range(height)]
    for line1, line2 in zip(horizontal_lines, horizontal_lines[1:]):
        new_edges = [(v1, v2) for v1, v2 in zip(line1[:-1], reversed(line2[1:]))]
        edges.extend(new_edges)
    return edges

def create_line_architecture(n_qubits, backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(n_qubits)
    edges = connect_vertices_in_line(vertices)
    graph.add_edges(edges)
    name = dynamic_size_architecture_name(LINE, n_qubits)
    return Architecture(name=name, coupling_graph=graph, backend=backend, **kwargs)

def create_circle_architecture(n_qubits, backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(n_qubits)
    edges = connect_vertices_in_line(vertices)
    edges.append((vertices[-1], vertices[0]))
    graph.add_edges(edges)
    name = dynamic_size_architecture_name(CIRCLE, n_qubits)
    return Architecture(name=name, coupling_graph=graph, backend=backend, **kwargs)

def create_square_architecture(n_qubits, backend=None, **kwargs):
    # No floating point errors
    sqrt_qubits = 0
    for n in range(n_qubits):
        if n_qubits == n**2:
            sqrt_qubits = n
        if n**2 > n_qubits:
            break
    if sqrt_qubits == 0:
        raise KeyError("Sqaure architecture requires a square number of qubits, but got " + str(n_qubits))
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(n_qubits)
    edges = connect_vertices_as_grid(sqrt_qubits, sqrt_qubits, vertices)
    graph.add_edges(edges)
    name = dynamic_size_architecture_name(SQUARE, n_qubits)
    return Architecture(name=name, coupling_graph=graph, backend=backend, **kwargs)

def create_ibm_qx2_architecture(**kwargs):
    m = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 1, 0, 1, 1],
        [0, 0, 1, 0, 1],
        [0, 0, 1, 1, 0]
    ])
    return Architecture(IBM_QX2, coupling_matrix=m, **kwargs)

def create_ibm_qx4_architecture(**kwargs):
    m = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 1, 0, 1, 1],
        [0, 0, 1, 0, 1],
        [0, 0, 1, 1, 0]
    ])
    return Architecture(IBM_QX4, coupling_matrix=m, **kwargs)

def create_ibm_qx3_architecture(**kwargs):
    m = np.array([
        #0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #0
        [1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #1
        [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #2
        [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #3
        [0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], #4
        [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], #5
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1], #6
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0], #7
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0], #8
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0], #9
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0], #10
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0], #11
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], #12
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0], #13
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1], #14
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]  #15
    ])
    return Architecture(IBM_QX3, coupling_matrix=m, **kwargs)

def create_ibm_qx5_architecture(**kwargs):
    m = np.array([
        #0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], #0
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], #1
        [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], #2
        [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], #3
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], #4
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], #5
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0], #6
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0], #7
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0], #8
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0], #9
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0], #10
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0], #11
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], #12
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0], #13
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1], #14
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]  #15
    ])
    return Architecture(IBM_QX5, coupling_matrix=m, **kwargs)

def create_ibm_q20_tokyo_architecture(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(20)
    edges = connect_vertices_as_grid(5, 4, vertices)
    cross_edges = [
        (1, 7), (2, 8),
        (3, 5), (4, 6),
        (6, 12), (7, 13),
        (8, 10), (9, 11),
        (11, 17), (12, 18),
        (13, 15), (14, 16)
    ]
    edges.extend([(vertices[v1], vertices[v2]) for v1, v2 in cross_edges])
    graph.add_edges(edges)
    return Architecture(name=IBM_Q20_TOKYO, coupling_graph=graph, backend=backend, **kwargs)

def create_ibmq_poughkeepsie(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(20)
    edges = connect_vertices_in_line(vertices)
    cross_edges = [
        (0,5), (7,1), (5,14), (10,19)
    ]
    edges.extend([(vertices[v1], vertices[v2]) for v1, v2 in cross_edges])
    graph.add_edges(edges)
    return Architecture(name=IBMQ_POUGHKEEPSIE, coupling_graph=graph, backend=backend, **kwargs)

def create_ibmq_singapore(backend=None, name=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(20)
    edges = connect_vertices_in_line([v for i, v in enumerate(vertices) if i not in [3, 15]])
    cross_edges = [(1,11), (3, 4), (5, 10), (9,14), (15, 16), (8,18)]
    edges.extend([(vertices[v1], vertices[v2]) for v1, v2 in cross_edges])
    graph.add_edges(edges)
    reduce_order = [19, 18, 17, 15, 16, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 3, 4, 2, 1, 0]
    if name is not None:
        return Architecture(name=IBMQ_BOEBLINGEN, coupling_graph=graph, backend=backend, reduce_order=reduce_order, **kwargs)
    return Architecture(name=IBMQ_SINGAPORE, coupling_graph=graph, backend=backend, reduce_order=reduce_order, **kwargs)

def create_rigetti_19q_acorn_architecture(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(20)
    edges = connect_vertices_in_line([vertices[i] for i in range(20) if i != 8])
    extra_edges = [(8,9), (0,18), (2,16), (4,14), (6, 12)]
    edges += [(vertices[v1], vertices[v2]) for v1, v2 in extra_edges]
    graph.add_edges(edges)
    return Architecture(RIGETTI_19Q_ACORN, coupling_graph=graph, reduce_order=[19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 8, 9, 7, 6, 5, 4, 3, 2, 1, 0], backend=backend, **kwargs)


def create_rigetti_16q_aspen_architecture(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(16)
    edges = connect_vertices_in_line(vertices)
    extra_edges = [(0, 7), (8, 15), (15, 0)]
    edges += [(vertices[v1], vertices[v2]) for v1, v2 in extra_edges]
    graph.add_edges(edges)
    return Architecture(RIGETTI_16Q_ASPEN, coupling_graph=graph, backend=backend, **kwargs)

def create_sycamore_like(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = list(graph.add_vertices(20))
    line = vertices[:1]+vertices[2:4]+vertices[5:13]+vertices[14:17]+vertices[18:]
    edges = connect_vertices_in_line(line)
    extra_edges = [(1,2),(2,6), (4,5), (4,14),(5,12), (6,11), (7,10), (13,14), (12,16), (11,18), (10,17), (17,18)]
    edges += [(vertices[v1], vertices[v2]) for v1, v2 in extra_edges]
    graph.add_edges(edges)
    return Architecture(SYCAMORE_LIKE, coupling_graph=graph, backend=backend, reduce_order=[19, 17, 18, 16,15,13,14,12,11,10,9,8,7,6,4,5,3,1,2,0], **kwargs)

def create_rigetti_8q_agave_architecture(**kwargs):
    m = np.array([
        [0, 1, 0, 0, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [1, 0, 0, 0, 0, 0, 1, 0]
    ])
    return Architecture(RIGETTI_8Q_AGAVE, coupling_matrix=m, **kwargs)

def create_recursive_architecture(**kwargs):
    m = np.array([
        #0  1  2  3  4  5  6  7  8
        [0, 0, 1, 0, 0, 0, 0, 0, 0],#0
        [0, 0, 1, 0, 0, 0, 0, 0, 0],#1
        [1, 1, 0, 0, 0, 0, 0, 1, 0],#2
        [0, 0, 0, 0, 0, 1, 0, 0, 0],#3
        [0, 0, 0, 0, 0, 1, 0, 0, 0],#4
        [0, 0, 0, 1, 1, 0, 0, 1, 0],#5
        [0, 0, 0, 0, 0, 0, 0, 1, 0],#6
        [0, 0, 1, 0, 0, 1, 1, 0, 1],#7
        [0, 0, 0, 0, 0, 0, 0, 1, 0] #8
    ])
    return Architecture(name=REC_ARCH, coupling_matrix=m, reduce_order=[8,6,4,3,5,7,1,2,0], **kwargs)

def create_fully_connected_architecture(n_qubits=None, **kwargs):
    if n_qubits is None:
        print("Warning: size is not given for the fully connected architecuture, using 9 as default.")
        n_qubits = 9
    m = np.ones(shape=(n_qubits, n_qubits))
    for i in range(n_qubits):
        m[i][i] = 0
    name = dynamic_size_architecture_name(FULLY_CONNECTED, n_qubits)
    return Architecture(name, coupling_matrix=m, **kwargs)

def create_dynamic_density_hamiltonian_architecture(n_qubits, density_prob=0.1, backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(n_qubits)
    edges = connect_vertices_in_line(vertices)
    n_edges = int(density_prob*n_qubits*(n_qubits-1)/2) - n_qubits+1 # Number of edges still to add.
    if n_edges > 0:
        possible_edges = [(v1, v2) for i, v1 in enumerate(vertices) for v2 in vertices[i+2:]]
        indices = np.random.choice(len(possible_edges), n_edges, replace=False)
        edges.extend([possible_edges[i] for i in indices])
    graph.add_edges(edges)
    name = dynamic_size_architecture_name(DENSITY+str(density_prob), n_qubits)
    return Architecture(name=name, coupling_graph=graph, backend=backend, **kwargs)

def create_dynamic_density_tree_architecture(n_qubits, density_prob=0.1, backend=None, **kwargs):
    graph = nx.Graph()
    vertices = [i for i in range(n_qubits)]
    graph.add_nodes_from(vertices)
    edges = []
    indices = [i for i in range(len(vertices))]
    index = np.random.choice(indices)
    root = vertices[index]
    indices.remove(index)
    stack = [root]
    while stack != []:
        parent = stack.pop(0)
        if indices != []:
            if len(indices) == 1:
                n_children = 1
            else:
                p = [0]
                x = 0.5
                while len(p) < len(indices)-1:
                    p.append(x)
                    x = x/2
                p.append(x*2)
                n_children = np.random.choice(len(indices), p=p) # Ensure that the parent has at least 1 child
            child_indices = np.random.choice(indices, n_children, replace=False)
            children = [vertices[child] for child in child_indices]
            edges += [(parent, child) for child in children]
            [indices.remove(i) for i in child_indices]
            stack += children
    n_edges = int(density_prob*n_qubits*(n_qubits-1)/2) - len(edges) # Number of edges still to add.
    if n_edges > 0:
        possible_edges = [(v1, v2) for i, v1 in enumerate(vertices) for v2 in vertices[i+1:] if (v1,v2) not in edges and v1!=v2 and (v2,v1) not in edges]
        indices = np.random.choice(len(possible_edges), n_edges, replace=False)
        edges.extend([possible_edges[i] for i in indices])
    #print(*sorted([(e1,e2) if e1 > e2 else (e2,e1) for e1,e2 in edges], key=lambda p: p[0]), sep="\n")
    graph.add_edges_from(edges)

    # Number the edges
    numbered = {}
    numbered = {v:i for i,v in enumerate(nx.dfs_postorder_nodes(graph, source=root))}
    # Make the coupling graph and adjust the numbering
    m = nx.to_numpy_array(graph)
    perm = [numbered[i] for i in range(n_qubits)]
    perm = [perm.index(i) for i in range(n_qubits)]
    m = m[perm][:, perm]
    spanning_tree = nx.dfs_tree(graph, source=root)
    #print(*[[numbered[n] for n in edge] for edge in spanning_tree.edges()])
    g = Graph(backend=backend)
    vertices = g.add_vertices(n_qubits)
    g.add_edges([(vertices[v1], vertices[v2]) for v1, v2 in [[numbered[n] for n in edge] for edge in spanning_tree.edges()]])
    arch = Architecture(name="temp", coupling_graph=g, backend=backend, **kwargs)
    subgraph = [i for i in range(n_qubits)]
    reduce_order = []
    while subgraph != []:
        non_cutting = arch.non_cutting_vertices(subgraph)
        node = max(non_cutting)
        reduce_order.append(node)
        subgraph.remove(node)
    #arch.reduce_order = reduce_order
    #reduce_order = []
    name = dynamic_size_architecture_name(DENSITY+str(density_prob), n_qubits)
    arch = Architecture(name, coupling_matrix=m, reduce_order=reduce_order, **kwargs)
    return arch

def create_dynamic_density_architecture(n_qubits, density_prob=0.1, backend=None, **kwargs):
    # Generate a random graph
    graph = nx.gnp_random_graph(n_qubits, density_prob)
    # Make sure it is connected
    connectivity = nx.all_pairs_node_connectivity(graph)
    disconnected = [(v1, v2) for v1,d in connectivity.items() for v2, score in d.items() if score == 0]
    while disconnected != []:
        to_add = np.random.choice(len(disconnected))
        graph.add_edge(*disconnected[to_add])
        connectivity = nx.all_pairs_node_connectivity(graph)
        disconnected = [(v1, v2) for v1,d in connectivity.items() for v2, score in d.items() if score == 0]
    # Number the edges
    numbered = {}
    numbered = {v:i for i,v in enumerate(nx.dfs_postorder_nodes(graph))}
    # Make the coupling graph and adjust the numbering
    m = nx.to_numpy_array(graph)
    perm = [numbered[i] for i in range(n_qubits)]
    perm = [perm.index(i) for i in range(n_qubits)]
    m = m[perm][:, perm]
    # Find the reduce order
    spanning_tree = nx.dfs_tree(graph)
    #print(*[[numbered[n] for n in edge] for edge in spanning_tree.edges()])
    g = Graph(backend=backend)
    vertices = g.add_vertices(n_qubits)
    g.add_edges([(vertices[v1], vertices[v2]) for v1, v2 in [[numbered[n] for n in edge] for edge in spanning_tree.edges()]])
    arch = Architecture(name="temp", coupling_graph=g, backend=backend, **kwargs)
    subgraph = [i for i in range(n_qubits)]
    reduce_order = []
    while subgraph != []:
        non_cutting = arch.non_cutting_vertices(subgraph)
        node = max(non_cutting)
        reduce_order.append(node)
        subgraph.remove(node)
    #arch.reduce_order = reduce_order
    #reduce_order = []
    name = dynamic_size_architecture_name(DENSITY+str(density_prob), n_qubits)
    arch = Architecture(name, coupling_matrix=m, reduce_order=reduce_order, **kwargs)
    return arch

def create_ibm_rochester(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(53)
    edges = connect_vertices_in_line([vertices[i] for i in range(53) if i not in [7, 14, 17,30,37, 40]])
    extra_edges = [(7, 8), (7,20), (14,15), (14,44), (17,18), (17, 28), (0,22), (30,31), (30,42), (37, 38), (40,41), (40,51)]
    edges += [(vertices[v1], vertices[v2]) for v1, v2 in extra_edges]
    graph.add_edges(edges)
    reduce_order = []
    i = 52
    for j in reversed([7, 14, 17,30,37, 40]):
        reduce_order += list(range(i, j+1, -1))
        reduce_order += [j, j+1]
        i = j-1
    reduce_order += list(range(i,-1, -1))
    #print(reduce_order)
    return Architecture(IBM_ROCHESTER, coupling_graph=graph, reduce_order=reduce_order, backend=backend, **kwargs)

def create_google_sycamore(backend=None, **kwargs):
    graph = Graph(backend=backend)
    vertices = graph.add_vertices(53)
    edges = connect_vertices_in_line([vertices[i] for i in range(53) if i not in [6,12,31,32,47,50]])
    extra_edges = []
    extra_edges += list(zip([0,1,4],[10,9,8]))
    extra_edges += list(zip(range(10, 5, -1),range(14,19)))
    extra_edges += list(zip(range(12,19),list(range(29,23, -1)) +[21]))
    extra_edges += list(zip(list(range(29,24, -1)), list(range(34,39))))
    extra_edges += list(zip(list(range(33,38)), [32] + list(range(43,39, -1))))
    extra_edges += list(zip(list(range(42,38, -1)), [45,46,48,47]))
    extra_edges += list(zip([45,46],[50,51]))
    extra_edges += list(zip([1,6,21,33,47,50,12,32,31],[4,7,24,38, 48,51,13,43,33]))
    edges += [(vertices[v1], vertices[v2]) for v1, v2 in extra_edges]
    graph.add_edges(edges)
    reduce_order = []
    i = 52
    for j in reversed([6,12,32, 47,50]):
        if j == 32:
            reduce_order += list(range(i, j+1, -1))
            reduce_order += [j,31, j+1]
            i = j-2
        else:
            reduce_order += list(range(i, j+1, -1))
            reduce_order += [j, j+1]
            i = j-1
    reduce_order += list(range(i,-1, -1))
    arch = Architecture(GOOGLE_SYCAMORE, coupling_graph=graph, reduce_order=reduce_order, backend=backend, **kwargs)
    return arch

def create_architecture(name, **kwargs):
    # Source Rigetti architectures: https://www.rigetti.com/qpu # TODO create the architectures from names in pyquil.list_quantum_computers() <- needs mapping
    # Source IBM architectures: http://iic.jku.at/files/eda/2018_tcad_mapping_quantum_circuit_to_ibm_qx.pdf​
    # IBM architectures are currently ignoring CNOT direction.
    if isinstance(name, Architecture):
        return name
    arch_dict = {}
    arch_dict[SYCAMORE_LIKE] = create_sycamore_like
    arch_dict[IBMQ_POUGHKEEPSIE] = create_ibmq_poughkeepsie
    arch_dict[IBMQ_SINGAPORE] = create_ibmq_singapore
    arch_dict[IBMQ_BOEBLINGEN] = lambda **kwargs : create_ibmq_singapore(name=IBMQ_BOEBLINGEN, **kwargs)
    arch_dict[DENSITY] = create_dynamic_density_tree_architecture
    arch_dict[GOOGLE_SYCAMORE] = create_google_sycamore
    arch_dict[IBM_ROCHESTER] = create_ibm_rochester
    if name == SQUARE:
        return create_square_architecture(**kwargs)
    elif name == LINE:
        return create_line_architecture(**kwargs)
    elif name == FULLY_CONNECTED:
        return create_fully_connected_architecture(**kwargs)
    elif name == CIRCLE:
        return create_circle_architecture(**kwargs)
    elif name == IBM_QX2:
        return create_ibm_qx2_architecture(**kwargs)
    elif name == IBM_QX3:
        return create_ibm_qx3_architecture(**kwargs)
    elif name == IBM_QX4:
        return create_ibm_qx4_architecture(**kwargs)
    elif name == IBM_QX5:
        return create_ibm_qx5_architecture(**kwargs)
    elif name == IBM_Q20_TOKYO:
        return create_ibm_q20_tokyo_architecture(**kwargs)
    elif name == RIGETTI_19Q_ACORN:
        return create_rigetti_19q_acorn_architecture(**kwargs)
    elif name == RIGETTI_16Q_ASPEN:
        return create_rigetti_16q_aspen_architecture(**kwargs)
    elif name == RIGETTI_8Q_AGAVE:
        return create_rigetti_8q_agave_architecture(**kwargs)
    elif name == REC_ARCH:
        return create_recursive_architecture(**kwargs)
    elif name in arch_dict.keys():
        return arch_dict[name](**kwargs)
    else:
        raise KeyError("name" + str(name) + "not recognized as architecture name. Please use one of", *architectures)

def colored_print_9X9(np_array):
    """
    Prints a 9x9 numpy array with colors representing their distance in a 9x9 square architecture
    :param np_array:  the array
    """
    if np_array.shape == (9,9):
        CRED = '\033[91m '
        CEND = '\033[0m '
        CGREEN = '\33[32m '
        CYELLOW = '\33[33m '
        CBLUE = '\33[34m '
        CWHITE = '\33[37m '
        CVIOLET = '\33[35m '
        color = [CBLUE, CGREEN, CVIOLET, CYELLOW, CRED]
        layout = [[0,1,2,3,2,1,2,3,4],
                  [1,0,1,2,1,2,3,2,3],
                  [2,1,0,1,2,3,4,3,2],
                  [3,2,1,0,1,2,3,2,1],
                  [2,1,2,1,0,1,2,1,2],
                  [1,2,3,2,1,0,1,2,3],
                  [2,3,4,3,2,1,0,1,2],
                  [3,2,3,2,1,2,1,0,1],
                  [4,3,2,1,2,3,2,1,0]]
        for i, l in enumerate(layout):
            print('[', ', '.join([(color[c] + '1' if v ==1 else CWHITE + '0') for c, v in zip(l, np_array[i])]), CEND, ']')
    else:
        print(np_array)

if __name__ == '__main__':
    sys.path.append('..')
    #n_qubits = 25
    #for name in dynamic_size_architectures:
    #    arch = create_architecture(name, n_qubits=n_qubits)
    #    arch.visualize()

    arch = create_architecture(RIGETTI_19Q_ACORN)
    arch.visualize()