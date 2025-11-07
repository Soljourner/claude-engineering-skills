# NetworkX Flow Networks - Reference Guide

## Graph Algorithms for Hydraulics

### 1. Path Finding Algorithms

#### Shortest Path Algorithms
Used to find minimum resistance paths, optimal routing, and critical flow paths.

**Simple Paths**
```python
# All simple paths between two nodes
paths = nx.all_simple_paths(G, source, target)

# Shortest path (unweighted)
path = nx.shortest_path(G, source, target)

# Shortest path (weighted by pipe length, resistance, etc.)
path = nx.shortest_path(G, source, target, weight='length')
length = nx.shortest_path_length(G, source, target, weight='length')

# All shortest paths (when multiple paths have same length)
paths = nx.all_shortest_paths(G, source, target, weight='resistance')
```

**Hydraulic Applications:**
- Find minimum pipe length routes
- Identify lowest resistance paths
- Optimize delivery routes in distribution networks
- Calculate friction head loss along paths

#### Dijkstra's Algorithm
```python
# Single source shortest paths to all nodes
lengths = nx.single_source_dijkstra_path_length(G, source, weight='resistance')
paths = nx.single_source_dijkstra_path(G, source, weight='resistance')

# Shortest path between specific nodes
length, path = nx.single_source_dijkstra(G, source, target, weight='resistance')
```

**Hydraulic Applications:**
- Minimize total pressure drop
- Find least-cost pipe routes
- Optimize pump placement

### 2. Flow Algorithms

#### Maximum Flow (Ford-Fulkerson / Edmonds-Karp)
```python
# Maximum flow from source to sink
flow_value, flow_dict = nx.maximum_flow(G, source, sink, capacity='capacity')

# Flow value only
flow_value = nx.maximum_flow_value(G, source, sink, capacity='capacity')

# Minimum cut
cut_value, partition = nx.minimum_cut(G, source, sink, capacity='capacity')
```

**Hydraulic Applications:**
- Determine maximum system capacity
- Identify bottleneck pipes (min cut)
- Size pipes for required flow
- Analyze redundancy and reliability
- Design expansion priorities

**Capacity Definition:**
For hydraulic networks, capacity can be defined as:
```python
# Based on velocity limit
capacity = A * v_max  # A = π(D/2)², v_max = 2-3 m/s

# Based on Reynolds number limit
capacity = (Re_max * μ) / (ρ * D)

# Based on pressure drop limit
capacity = sqrt((2 * D * ΔP_max) / (ρ * f * L))
```

#### Minimum Cost Flow
```python
# Find flow that satisfies demands at minimum cost
flowDict = nx.min_cost_flow(G)
flowCost = nx.cost_of_flow(G, flowDict)

# With node demands
G.nodes['J1']['demand'] = -100  # Source (negative)
G.nodes['J2']['demand'] = 50    # Sink (positive)
G.nodes['J3']['demand'] = 50    # Sink (positive)
```

**Hydraulic Applications:**
- Minimize pumping costs
- Optimize energy consumption
- Balance flow distribution
- Pump scheduling in multi-source systems

**Cost Definition:**
```python
# Energy cost (proportional to head loss)
cost = f * L / D  # Friction factor × length / diameter

# Pumping energy cost
cost = ρ * g * h * Q * time / η  # Where h = head, η = efficiency

# Pipe material cost
cost = price_per_kg * (π/4 * (D_outer² - D_inner²)) * L * ρ_material
```

### 3. Connectivity and Reliability

#### Network Connectivity
```python
# Check if network is connected
is_connected = nx.is_weakly_connected(G)  # For directed graphs
is_connected = nx.is_connected(G.to_undirected())

# Node connectivity (min nodes to remove to disconnect)
k = nx.node_connectivity(G, source, target)

# Edge connectivity (min edges to remove to disconnect)
k = nx.edge_connectivity(G, source, target)

# Find articulation points (critical junctions)
articulation_points = list(nx.articulation_points(G.to_undirected()))

# Find bridges (critical pipes)
bridges = list(nx.bridges(G.to_undirected()))
```

**Hydraulic Applications:**
- Identify single points of failure
- Assess network redundancy
- Plan maintenance shutdowns
- Design backup flow paths
- Emergency response planning

#### Component Analysis
```python
# Find connected components
components = list(nx.connected_components(G.to_undirected()))
n_components = nx.number_connected_components(G.to_undirected())

# Strongly connected components (for directed graphs)
strong_components = list(nx.strongly_connected_components(G))
```

**Hydraulic Applications:**
- Identify isolated subsystems
- Verify network completeness
- Detect disconnected zones

### 4. Network Topology Analysis

#### Degree Analysis
```python
# Node degree (number of connections)
degree_dict = dict(G.degree())

# In-degree and out-degree for directed graphs
in_degree = dict(G.in_degree())
out_degree = dict(G.out_degree())

# Find hub junctions
hubs = [n for n, d in G.degree() if d >= threshold]
```

**Hydraulic Applications:**
- Identify major junctions (high degree nodes)
- Find terminal points (degree = 1)
- Analyze network complexity

#### Centrality Measures
```python
# Betweenness centrality (how often node appears on shortest paths)
betweenness = nx.betweenness_centrality(G, weight='resistance')

# Closeness centrality (average distance to all other nodes)
closeness = nx.closeness_centrality(G, distance='length')

# Degree centrality (fraction of nodes connected to)
degree_centrality = nx.degree_centrality(G)
```

**Hydraulic Applications:**
- Identify critical junctions for monitoring
- Optimize sensor placement
- Prioritize maintenance and upgrades
- Assess vulnerability to failures

### 5. Cycle Detection

#### Finding Cycles
```python
# Find all simple cycles
cycles = list(nx.simple_cycles(G))

# Check if graph is acyclic
is_dag = nx.is_directed_acyclic_graph(G)

# Find cycle basis (for undirected graphs)
cycles = nx.cycle_basis(G.to_undirected())
```

**Hydraulic Applications:**
- Identify loop flows
- Validate network topology
- Analyze flow distribution in meshed networks
- Hardy Cross method applications

### 6. Graph Traversal

#### Depth-First Search (DFS)
```python
# DFS tree
dfs_tree = nx.dfs_tree(G, source)

# DFS edges
dfs_edges = list(nx.dfs_edges(G, source))

# Topological sort (for DAGs)
if nx.is_directed_acyclic_graph(G):
    topo_order = list(nx.topological_sort(G))
```

**Hydraulic Applications:**
- Trace flow paths systematically
- Order calculation for sequential systems
- Generate system schematic from connections

#### Breadth-First Search (BFS)
```python
# BFS tree
bfs_tree = nx.bfs_tree(G, source)

# BFS edges
bfs_edges = list(nx.bfs_edges(G, source))
```

**Hydraulic Applications:**
- Find nodes at specific distances
- Analyze propagation of pressure waves
- Zone identification

### 7. Subgraph Analysis

#### Finding Subgraphs
```python
# Subgraph induced by nodes
nodes_subset = ['J1', 'J2', 'J3']
subgraph = G.subgraph(nodes_subset)

# Edge-induced subgraph
edges_subset = [('J1', 'J2'), ('J2', 'J3')]
edge_subgraph = G.edge_subgraph(edges_subset)

# k-core (subgraph where all nodes have degree >= k)
k_core = nx.k_core(G, k=2)
```

**Hydraulic Applications:**
- Analyze subsystems independently
- Isolate zones for maintenance
- Study local network characteristics

## NetworkX Documentation

### Official Documentation
- **Main Documentation**: https://networkx.org/documentation/stable/
- **Tutorial**: https://networkx.org/documentation/stable/tutorial.html
- **Reference**: https://networkx.org/documentation/stable/reference/index.html

### Specific Algorithm References

#### Flow Algorithms
- **Maximum Flow**: https://networkx.org/documentation/stable/reference/algorithms/flow.html
- **Min Cost Flow**: https://networkx.org/documentation/stable/reference/algorithms/flow.html#minimum-cost-flow
- **Capacity Scaling**: https://networkx.org/documentation/stable/reference/algorithms/flow.html#capacity-scaling-minimum-cost-flow

#### Path Algorithms
- **Shortest Paths**: https://networkx.org/documentation/stable/reference/algorithms/shortest_paths.html
- **All Pairs Shortest Path**: https://networkx.org/documentation/stable/reference/algorithms/shortest_paths.html#all-pairs-shortest-paths
- **Dijkstra**: https://networkx.org/documentation/stable/reference/algorithms/shortest_paths.html#dijkstra-s-algorithm

#### Connectivity
- **Connectivity**: https://networkx.org/documentation/stable/reference/algorithms/connectivity.html
- **Components**: https://networkx.org/documentation/stable/reference/algorithms/component.html
- **Cuts**: https://networkx.org/documentation/stable/reference/algorithms/cuts.html

#### Graph Theory
- **Cycles**: https://networkx.org/documentation/stable/reference/algorithms/cycles.html
- **Centrality**: https://networkx.org/documentation/stable/reference/algorithms/centrality.html
- **Traversal**: https://networkx.org/documentation/stable/reference/algorithms/traversal.html

### Data Structures
- **Graph**: https://networkx.org/documentation/stable/reference/classes/graph.html
- **DiGraph**: https://networkx.org/documentation/stable/reference/classes/digraph.html
- **MultiGraph**: https://networkx.org/documentation/stable/reference/classes/multigraph.html
- **MultiDiGraph**: https://networkx.org/documentation/stable/reference/classes/multidigraph.html

### Import/Export
- **GraphML**: https://networkx.org/documentation/stable/reference/readwrite/graphml.html
- **GML**: https://networkx.org/documentation/stable/reference/readwrite/gml.html
- **JSON**: https://networkx.org/documentation/stable/reference/readwrite/json_graph.html
- **Adjacency List**: https://networkx.org/documentation/stable/reference/readwrite/adjlist.html

## Related Resources

### Graph Theory for Engineering
- **Graph Theory and Applications** (Gross & Yellen)
- **Networks, Crowds, and Markets** (Easley & Kleinberg) - Free online: https://www.cs.cornell.edu/home/kleinber/networks-book/
- **Introduction to Graph Theory** (West)

### Hydraulic Network Analysis
- **Water Distribution Systems Handbook** (Mays)
- **Fundamentals of Hydraulic Engineering Systems** (Hwang & Houghtalen)
- **Applied Hydraulics in Engineering** (Roberson & Crowe)

### Python Scientific Computing
- **NumPy Documentation**: https://numpy.org/doc/stable/
- **SciPy Documentation**: https://docs.scipy.org/doc/scipy/
- **Matplotlib Documentation**: https://matplotlib.org/stable/contents.html

### Specialized Libraries

#### EPANET (Hydraulic Network Solver)
- **EPANET**: https://www.epa.gov/water-research/epanet
- **wntr (Water Network Tool for Resilience)**: https://wntr.readthedocs.io/
  - Python wrapper for EPANET
  - NetworkX integration for network analysis
  - Resilience metrics and analysis

#### Graph Visualization
- **Graphviz**: https://graphviz.org/
- **PyGraphviz**: https://pygraphviz.github.io/
- **Plotly Network Graphs**: https://plotly.com/python/network-graphs/

## Quick Reference Table

| Task | NetworkX Function | Hydraulic Application |
|------|------------------|----------------------|
| Max flow capacity | `nx.maximum_flow()` | System capacity analysis |
| Bottleneck identification | `nx.minimum_cut()` | Find limiting pipes |
| Shortest path | `nx.shortest_path()` | Minimum resistance route |
| Critical pipes | `nx.bridges()` | Single point of failure |
| Critical junctions | `nx.articulation_points()` | Single point of failure |
| All flow paths | `nx.all_simple_paths()` | Parallel flow analysis |
| Node importance | `nx.betweenness_centrality()` | Sensor placement |
| Network loops | `nx.simple_cycles()` | Hardy Cross analysis |
| Connectivity check | `nx.is_connected()` | System validation |
| Degree analysis | `G.degree()` | Junction complexity |

## Common Patterns

### Pattern 1: Hydraulic Resistance Network
```python
# Set resistance as edge weight
for u, v, data in G.edges(data=True):
    f = 0.02  # friction factor
    data['resistance'] = f * data['L'] / data['D']

# Find minimum resistance path
path = nx.shortest_path(G, source, sink, weight='resistance')
total_resistance = nx.shortest_path_length(G, source, sink, weight='resistance')
```

### Pattern 2: Multi-Commodity Flow
```python
# Multiple sources and sinks
sources = {'Tank1': 50, 'Tank2': 30}  # Supply
sinks = {'Outlet1': 20, 'Outlet2': 35, 'Outlet3': 25}  # Demand

# Solve for each commodity
for source, supply in sources.items():
    for sink, demand in sinks.items():
        flow = nx.maximum_flow_value(G, source, sink)
        allocation = min(supply, demand, flow)
```

### Pattern 3: Network Reliability
```python
# Calculate network reliability (all paths must exist)
def network_reliability(G, source, sink):
    try:
        # Check connectivity
        if not nx.has_path(G, source, sink):
            return 0.0

        # Count independent paths
        edge_connectivity = nx.edge_connectivity(G, source, sink)

        # Reliability increases with redundancy
        reliability = 1 - (1 / (1 + edge_connectivity))
        return reliability
    except:
        return 0.0
```

### Pattern 4: Pressure Drop Calculation
```python
# Calculate pressure drop along path
def calculate_pressure_drop(G, path, Q, rho=1000, mu=0.001):
    """Calculate total pressure drop along flow path"""
    total_dP = 0

    for u, v in zip(path[:-1], path[1:]):
        edge = G[u][v]
        D, L = edge['D'], edge['L']
        A = np.pi * (D/2)**2
        v = Q / A  # velocity

        # Reynolds number
        Re = rho * v * D / mu

        # Friction factor (Colebrook-White approximation)
        if Re > 2300:  # Turbulent
            f = 0.02  # Simplified
        else:  # Laminar
            f = 64 / Re

        # Darcy-Weisbach equation
        dP = f * (L/D) * (rho * v**2 / 2)
        total_dP += dP
        edge['dP'] = dP

    return total_dP
```

## Tips and Best Practices

1. **Choose the right graph type:**
   - Use `DiGraph` for flow networks (flow has direction)
   - Use `MultiDiGraph` for parallel pipes between same junctions
   - Use `Graph` for undirected connectivity analysis

2. **Consistent units:** Always use SI units (m, m³/s, Pa, kg/m³)

3. **Edge attributes:** Store all relevant pipe properties (D, L, k, material, age)

4. **Node attributes:** Store junction properties (elevation, demand, pressure)

5. **Validation:** Always check connectivity before flow analysis

6. **Performance:** For large networks (>1000 nodes), use sparse algorithms

7. **Visualization:** Use hierarchical layouts for tree-like networks

8. **Documentation:** Comment edge attribute meanings and units

## See Also

- `SKILL.md` - Complete guide with examples
- `examples.py` - Working code examples
- EPANET for detailed hydraulic simulation
- SciPy for numerical solvers
