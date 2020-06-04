from collections import deque


def flightLegs(Alist,start,dest):
    """
    Find the minimum number of flights required to travel between start and dest,
    and  determine the number of distinct routes between start and dest which
    require this minimum number of flights.
    Input:
        Alist: Adjacency list for airport network
        start, dest: starting and destination airports for journey.
        Airports are numbered from 0 to N-1 (inclusive) where N = len(Alist)
    Output:
        Flights: 2-element list containing:
        [the min number of flights between start and dest, the number of distinct
        jouneys with this min number]
        Return an empty list if no journey exist which connect start and dest
    """

    # Initialization
    L2 = [0 for l in Alist]  # Labels
    L3 = L2.copy()  # Distances when first discovered (also the minimum distance with bfs)
    L4 = L2.copy()  # Number of shortest paths
    Q = deque([start])
    L2[start] = 1
    L4[start] = 1

    # Loop thorugh all connected nodes to start node
    while Q:
        x = Q.popleft() # Remove node from front of queue

        for v in Alist[x]:
            if L2[v] == 0:
                # Add unexplored neighbors to back of queue
                Q.append(v)
                L2[v] = 1
                L3[v] = 1+L3[x]

            # Check whether this is a shortest path to v
            if L3[x]+1 == L3[v]:
                # Add number of shortest paths to x to number of shortest paths to v
                L4[v] += L4[x]

        # Check whether dest was found during iteration of x
        if L2[dest]:
            depth = L3[dest]
            # Iterate through all neighbours of dest and terminate
            for v in Alist[dest]:
                if L3[v]+1 == depth:
                    L4[dest] += L4[v]

            # Subtract duplicate count from node x, which is a neighbour of dest!!!
            L4[dest] -= L4[x]
            return [L3[dest], L4[dest]]

    # dest node not connected to start node
    return []


def safeJourney(Alist,start,dest):
    """
    Find safest journey from station start to dest
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the density for the connection.
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing safest journey and safety factor for safest journey
    """

    # Initialize dictionaries
    dinit = float('inf')
    Edict = {} # Explored nodes
    Udict = {} # Unexplored nodes
    N = len(Alist)

    Udict = {i: [None, dinit] for i in range(N)}
    Udict[start][1] = 0

    # Main search
    while Udict:
        # Find node with min d in Udict and move to Edict
        dmin = dinit
        nmin = None
        for n, v in Udict.items():
            d = v[1]
            if d < dmin:
                dmin = d
                nmin = n

        if nmin is None:
            # Destination node is disconnected
            return []

        # Add node to explored
        Edict[nmin] = Udict.pop(nmin)

        if nmin == dest:
            # Solution is complete
            # Route from start to dest
            path = deque([dest])
            node = dest
            while node != start:
                node = Edict[node][0]
                path.appendleft(node)
            return [list(path), Edict[dest][1]]

        # Update provisional safeties for unexplored neighbors of nmin
        for n, d in Alist[nmin]:
            if n in Udict:
                dcomp = max(d, dmin)
                if dcomp < Udict[n][1]:
                    # Replace distance and last node with better one
                    Udict[n] = [nmin, dcomp]

    return []


def shortJourney(Alist,start,dest):
    """
    Find shortest journey from station start to dest. If multiple shortest journeys
    exist, select journey which goes through the smallest number of stations.
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the time for the connection (rounded to the nearest minute).
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing shortest journey and duration of shortest journey
    """

    # Initialize dictionaries
    dinit = float('inf')
    Edict = {} # Explored nodes
    Udict = {} # Unexplored nodes
    N = len(Alist)

    Udict = {i: [None, dinit, 1] for i in range(N)}
    Udict[start][1] = 0

    # Main search
    while Udict:
        # Find node with min d in Udict and move to Edict
        dmin = dinit
        nmin = None
        for n, v in Udict.items():
            d = v[1]
            if d < dmin:
                dmin = d
                nmin = n

        if nmin is None:
            # Destination node is disconnected
            return []

        # Add node to explored
        Edict[nmin] = Udict.pop(nmin)

        if nmin == dest:
            # Solution is complete
            # Route from start to dest
            length = Edict[dest][2]
            path = [dest] * length
            node = dest
            for i in range(length-2, -1, -1):
                node = Edict[node][0]
                path[i] = node
            return [path, Edict[dest][1]]

        # Update times distances for unexplored neighbors of nmin
        for n, d in Alist[nmin]:
            if n in Udict:
                dcomp = dmin + d
                if dcomp < Udict[n][1]:
                    # Replace distance and last node with better one
                    Udict[n] = [nmin, dcomp, Edict[nmin][2]+1]
                elif dcomp == Udict[n][1]:
                    # Choose journey with smallest number of steps
                    new_length = Edict[nmin][2]+1
                    old_length = Udict[n][2]
                    if new_length < old_length:
                        Udict[n] = [nmin, dcomp, new_length]

    return []


def cheapCycling(SList,CList):
    """
    Find first and last stations for cheapest cycling trip
    Input:
        Slist: list whose ith element contains cheapest fare for arrival at and
        return from the ith station (stored in a 2-element list or tuple)
        Clist: list whose ith element contains a list of stations which can be
        cycled to directly from station i
    Stations are numbered from 0 to N-1 with N = len(Slist) = len(Clist)
    Output:
        stations: two-element list containing first and last stations of journey
    """

    N = len(CList)
    V = {i for i in range(N)}
    Q = set()
    total_cost = float('inf')
    nodes = [None, None]
    W = dict() # Connections
    counter = -1

    while V:
        # Take any unexplored node
        y = V.pop()
        if CList[y] == []:
            continue

        Q.add(y)
        counter += 1
        W[counter] = set([y])

        # Run until all nodes have been searched
        while Q:
            x = Q.pop() # Remove any node from set

            for v in CList[x]:
                if v not in W[counter]:
                    W[counter].add(v)
                    V.remove(v)
                    Q.add(v) # Add unexplored neighbors to set

    for k, conn_graph in W.items():
        admin = float('inf')
        admin2 = float('inf')
        anmin = None
        anmin2 = None
        rdmin = float('inf')
        rdmin2 = float('inf')
        rnmin = None
        rnmin2 = None
        # Calculate minimum arrival costs
        for node in conn_graph:
            d = SList[node][0]
            if d < admin2:
                if d < admin:
                    admin2 = admin
                    anmin2 = anmin
                    admin = d
                    anmin = node
                else:
                    admin2 = d
                    anmin2 = node

        # Calculate minimum return costs
        for node in conn_graph:
            d = SList[node][1]
            if d < rdmin2:
                if d < rdmin:
                    rdmin2 = rdmin
                    rnmin2 = rnmin
                    rdmin = d
                    rnmin = node
                else:
                    rdmin2 = d
                    rnmin2 = node

        # Compare cost to other connected graphs
        if anmin == rnmin:
            # Arrival and return nodes not distinct
            new_cost1 = admin2 + rdmin
            new_cost2 = admin + rdmin2

            # Code readability favoured over minute efficiency gains with extra if loops
            if new_cost1 < total_cost:
                total_cost = new_cost1
                nodes = [anmin2, rnmin]

            if new_cost2 < total_cost:
                total_cost = new_cost2
                nodes = [anmin, rnmin2]

        else:
            new_cost = admin + rdmin
            if new_cost < total_cost:
                total_cost = new_cost
                nodes = [anmin, rnmin]

    return nodes


if __name__=='__main__':
    pass
