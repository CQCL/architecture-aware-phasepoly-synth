from pyzx.linalg import Mat2

debug = False

def rec_steiner_gauss(matrix, architecture, full_reduce=False, x=None, y=None, permutation=None, **kwargs):
    """
    Performs recursive Gaussian elimination that is constraint bij the given architecture
    
    :param matrix: PyZX Mat2 matrix to be reduced
    :param architecture: The Architecture object to conform to
    :param full_reduce: Whether to fully reduce or only create an upper triangular form
    :param x: 
    :param y: 
    """
    #print(matrix)
    if permutation is None:
        permutation = [i for i in range(len(matrix.data))]
    else:
        matrix = Mat2([[row[i] for i in permutation] for row in matrix.data])
    #print(matrix)
    def row_add(c0, c1):
        matrix.row_add(c0, c1)
        debug and print("Reducing", c0, c1)
        c0 = architecture.qubit_map[c0]
        c1 = architecture.qubit_map[c1]
        if x != None: x.row_add(c0, c1)
        if y != None: y.col_add(c1, c0)
    def steiner_reduce(col, root, nodes, usable_nodes, rec_nodes, upper):
        generator = steiner_reduce_column(architecture, [row[col] for row in matrix.data], root, nodes, usable_nodes, rec_nodes, upper)
        cnot = next(generator, None)
        while cnot is not None:
            row_add(*cnot)
            cnot = next(generator, None)
    def rec_step(cols, rows):
        size = len(rows)
        # Upper triangular form uses the same structure.
        p_cols = []
        pivot = 0
        rows2 = [r for r in range(len(matrix.data[0])) if r in rows]
        cols2 = [c for c in range(len(matrix.data)) if c in cols]
        for i, c in enumerate(cols2):
            if pivot < size:
                nodes = [r for r in rows2[pivot:] if rows2[pivot]==r or matrix.data[r][c] == 1]
                steiner_reduce(c, rows2[pivot], nodes, cols2[i:], [], True)
                if matrix.data[rows2[pivot]][c] == 1:
                    p_cols.append(c)
                    pivot += 1
        # Full reduce requires the recursion
        if full_reduce:
            pivot -= 1
            for i, c in enumerate(cols):
                if c in p_cols:
                    nodes = [r for r in rows if r==rows[pivot] or matrix.data[r][c] == 1]
                    usable_nodes = cols[i:]
                    #rec_nodes = list(set([node for edge in architecture.distances["upper"][c][(max(usable_nodes),c)][1] for node in edge]))
                    #rec_nodes = [n for n in usable_nodes if n in rec_nodes]
                    rec_nodes = architecture.shortest_path(c, max(usable_nodes), usable_nodes) 

                    if len(nodes) > 1:
                        #print("-", c, nodes, cols, rec_nodes)
                        steiner_reduce(c, rows[pivot], nodes, cols, rec_nodes, False)
                        
                    # Do recursion on the given nodes.
                    if len(rec_nodes) > 1:
                        rec_step(list(reversed(rec_nodes)), rec_nodes)
                pivot -= 1
        #return rank
    qubit_order = architecture.reduce_order
    rec_step(qubit_order, list(reversed(qubit_order)))
    #print("Done!")

def steiner_reduce_column(architecture, col, root, nodes, usable_nodes, rec_nodes, upper):
    steiner_tree = architecture.rec_steiner_tree(root, nodes, usable_nodes, rec_nodes, upper)
    # Remove all zeros
    next_check = next(steiner_tree)
    debug and print("Step 1: remove zeros")
    if upper:
        zeros = []
        while next_check is not None:
            s0, s1 = next_check
            if col[s0] == 0:  # s1 is a new steiner point or root = 0
                zeros.append(next_check)
            next_check = next(steiner_tree)
        while len(zeros) > 0:
            # Iterate through zeros in reverse, because we know that all leaves are 1.
            s0, s1 = zeros.pop(-1)
            if col[s0] == 0:
                col[s0] = (col[s1]+col[s0])%2
                yield s1, s0
                debug and print(col[s0], col[s1])
    else:
        debug and print("deal with zero root")
        if next_check is not None and col[next_check[0]] == 0:  # root is zero
            print("WARNING : Root is 0 => reducing non-pivot column", col)
        debug and print("Step 1: remove zeros", col)
        while next_check is not None:
            s0, s1 = next_check
            if col[s1] == 0:  # s1 is a new steiner point
                col[s1] = (col[s1]+col[s0])%2
                yield s0, s1
            next_check = next(steiner_tree)
    # Reduce stuff
    debug and print("Step 2: remove ones")
    next_add = next(steiner_tree)
    while next_add is not None:
        s0, s1 = next_add
        col[s1] = (col[s1]+col[s0])%2
        yield s0, s1
        next_add = next(steiner_tree)
        debug and print(next_add)
    debug and print("Step 3: profit")