---
name: networkx-flow-networks
description: "Model hydraulic networks and multi-pump systems using graph theory"
category: packages
domain: fluids
complexity: intermediate
dependencies: [networkx, numpy]
---

# NetworkX Flow Networks

Master graph theory and network flow algorithms for modeling hydraulic piping systems, multi-pump configurations, and complex flow distribution networks.

## Installation

```bash
pip install networkx numpy matplotlib
```

## Graph Basics

NetworkX provides powerful tools for creating and analyzing directed/undirected graphs:

```python
import networkx as nx
import numpy as np

# Create a directed graph for flow networks
G = nx.DiGraph()

# Add nodes (junctions, pumps, tanks)
G.add_node('Tank', type='source', elevation=30)
G.add_node('J1', type='junction', elevation=0)
G.add_node('J2', type='junction', elevation=0)
G.add_node('Outlet', type='sink', elevation=0)

# Add edges (pipes) with attributes
G.add_edge('Tank', 'J1', diameter=0.15, length=100, roughness=0.00015)
G.add_edge('J1', 'J2', diameter=0.10, length=50, roughness=0.00015)
G.add_edge('J2', 'Outlet', diameter=0.10, length=75, roughness=0.00015)
```

## Engineering Applications

### Piping Network Representation

Model complex piping systems as directed graphs where:
- **Nodes** represent junctions, tanks, pumps, valves
- **Edges** represent pipes with flow capacity
- **Edge attributes** store diameter, length, roughness, resistance

```python
def create_piping_network():
    """Create a directed graph representing a piping network"""
    G = nx.DiGraph()

    # Add nodes with elevation data
    nodes = {
        'Reservoir': {'type': 'source', 'elevation': 50, 'pressure': 0},
        'Pump1': {'type': 'pump', 'elevation': 0},
        'J1': {'type': 'junction', 'elevation': 0},
        'J2': {'type': 'junction', 'elevation': 0},
        'Tank': {'type': 'sink', 'elevation': 20, 'demand': 100}
    }

    for node, attrs in nodes.items():
        G.add_node(node, **attrs)

    # Add pipes with hydraulic properties
    pipes = [
        ('Reservoir', 'Pump1', {'D': 0.2, 'L': 50, 'k': 0.00015}),
        ('Pump1', 'J1', {'D': 0.15, 'L': 100, 'k': 0.00015}),
        ('J1', 'J2', {'D': 0.10, 'L': 80, 'k': 0.00015}),
        ('J2', 'Tank', {'D': 0.10, 'L': 60, 'k': 0.00015})
    ]

    for u, v, attrs in pipes:
        G.add_edge(u, v, **attrs)

    return G
```

### Flow Distribution in Parallel Systems

Analyze flow distribution in networks with parallel paths:

```python
def analyze_parallel_flow(G, source, sink):
    """Analyze flow distribution in parallel pipe systems"""
    # Find all simple paths from source to sink
    all_paths = list(nx.all_simple_paths(G, source, sink))

    print(f"Number of flow paths: {len(all_paths)}")
    for i, path in enumerate(all_paths, 1):
        print(f"Path {i}: {' -> '.join(path)}")

        # Calculate total resistance for each path
        total_length = sum(G[u][v].get('L', 0) for u, v in zip(path[:-1], path[1:]))
        print(f"  Total length: {total_length} m")

    return all_paths

# Calculate resistance for Darcy-Weisbach
def pipe_resistance(D, L, f=0.02):
    """Calculate pipe resistance coefficient K = f*L/D"""
    return f * L / D

# Add resistance to edges
G = create_piping_network()
for u, v, data in G.edges(data=True):
    data['K'] = pipe_resistance(data['D'], data['L'])
```

### Path Analysis

Find critical paths, shortest paths, and identify bottlenecks:

```python
def path_analysis(G, source, sink):
    """Perform comprehensive path analysis"""

    # Shortest path by number of pipes
    shortest_path = nx.shortest_path(G, source, sink)
    print(f"Shortest path (by hops): {' -> '.join(shortest_path)}")

    # Shortest path by total length
    for u, v, data in G.edges(data=True):
        data['weight'] = data['L']

    shortest_length_path = nx.shortest_path(G, source, sink, weight='weight')
    total_length = nx.shortest_path_length(G, source, sink, weight='weight')
    print(f"Shortest path (by length): {' -> '.join(shortest_length_path)}")
    print(f"Total length: {total_length} m")

    # Find bottlenecks (minimum diameter in path)
    min_diameter = min(G[u][v]['D'] for u, v in zip(shortest_path[:-1], shortest_path[1:]))
    print(f"Bottleneck diameter: {min_diameter} m")

    return shortest_path, shortest_length_path
```

### Network Optimization

Optimize flow distribution and identify critical components:

```python
def find_critical_pipes(G):
    """Identify critical pipes whose removal would disconnect the network"""
    critical_edges = list(nx.bridges(G.to_undirected()))

    print("Critical pipes (bridges):")
    for u, v in critical_edges:
        data = G[u][v] if G.has_edge(u, v) else G[v][u]
        print(f"  {u} <-> {v}: D={data['D']} m, L={data['L']} m")

    return critical_edges

def redundancy_analysis(G, source, sink):
    """Analyze network redundancy"""
    # Node connectivity (minimum nodes to remove to disconnect)
    node_connectivity = nx.node_connectivity(G, source, sink)
    print(f"Node connectivity: {node_connectivity}")

    # Edge connectivity (minimum edges to remove to disconnect)
    edge_connectivity = nx.edge_connectivity(G, source, sink)
    print(f"Edge connectivity: {edge_connectivity}")

    if edge_connectivity > 1:
        print("Network has redundant paths")
    else:
        print("Network has no redundancy - single point of failure")
```

### Max Flow / Min Cut

Apply max flow algorithms for capacity analysis:

```python
def max_flow_analysis(G, source, sink):
    """Calculate maximum flow capacity using Ford-Fulkerson"""

    # Set capacity based on pipe diameter (Q ∝ D^2 for velocity limit)
    for u, v, data in G.edges(data=True):
        # Capacity proportional to cross-sectional area
        data['capacity'] = np.pi * (data['D']/2)**2 * 2.0  # Assuming 2 m/s max velocity

    # Calculate max flow
    flow_value, flow_dict = nx.maximum_flow(G, source, sink, capacity='capacity')

    print(f"Maximum flow capacity: {flow_value:.4f} m³/s")
    print("\nFlow distribution:")
    for u in flow_dict:
        for v, flow in flow_dict[u].items():
            if flow > 0:
                capacity = G[u][v]['capacity']
                utilization = flow / capacity * 100
                print(f"  {u} -> {v}: {flow:.4f} m³/s ({utilization:.1f}% capacity)")

    # Find minimum cut
    cut_value, partition = nx.minimum_cut(G, source, sink, capacity='capacity')
    reachable, non_reachable = partition

    print(f"\nMinimum cut value: {cut_value:.4f} m³/s")
    print(f"Cut set (bottleneck pipes):")
    for u, v, data in G.edges(data=True):
        if u in reachable and v in non_reachable:
            print(f"  {u} -> {v}: D={data['D']} m, capacity={data['capacity']:.4f} m³/s")

    return flow_value, flow_dict
```

## Hydraulic Network Modeling

Complete workflow for hydraulic network analysis:

```python
def hydraulic_network_model(G, source, sink, fluid_props):
    """Complete hydraulic network analysis"""

    print("=== Network Topology ===")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Pipes: {G.number_of_edges()}")
    print(f"Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")

    # Check connectivity
    if not nx.is_connected(G.to_undirected()):
        print("WARNING: Network is not fully connected!")
        components = list(nx.connected_components(G.to_undirected()))
        print(f"Number of separate components: {len(components)}")

    print("\n=== Path Analysis ===")
    paths = analyze_parallel_flow(G, source, sink)

    print("\n=== Critical Components ===")
    find_critical_pipes(G)

    print("\n=== Redundancy Analysis ===")
    redundancy_analysis(G, source, sink)

    print("\n=== Max Flow Analysis ===")
    max_flow_analysis(G, source, sink)

    # Calculate pressure drops (simplified)
    print("\n=== Pressure Drop Analysis ===")
    for u, v, data in G.edges(data=True):
        rho = fluid_props['density']
        mu = fluid_props['viscosity']
        f = 0.02  # Friction factor (simplified)

        # Estimate velocity based on diameter
        V = 2.0  # m/s assumed
        dP = f * (data['L'] / data['D']) * (rho * V**2 / 2)
        data['dP'] = dP

        print(f"{u} -> {v}: ΔP = {dP/1000:.2f} kPa")

# Example usage
G = create_piping_network()
fluid_props = {'density': 1000, 'viscosity': 0.001}  # Water at 20°C
hydraulic_network_model(G, 'Reservoir', 'Tank', fluid_props)
```

### Multi-Pump System Modeling

Model systems with multiple pumps in series or parallel:

```python
def create_multipump_network():
    """Create network with parallel and series pumps"""
    G = nx.DiGraph()

    # Source tank
    G.add_node('Source', type='tank', elevation=0)

    # Parallel pump configuration
    G.add_node('Pump1A', type='pump', head=50, flow_rated=0.05)
    G.add_node('Pump1B', type='pump', head=50, flow_rated=0.05)
    G.add_node('Junction1', type='junction')

    # Series pump configuration
    G.add_node('Pump2', type='pump', head=30, flow_rated=0.08)
    G.add_node('Junction2', type='junction')

    # Destination
    G.add_node('Destination', type='tank', elevation=60)

    # Connect parallel pumps
    G.add_edge('Source', 'Pump1A', D=0.15, L=10)
    G.add_edge('Source', 'Pump1B', D=0.15, L=10)
    G.add_edge('Pump1A', 'Junction1', D=0.15, L=20)
    G.add_edge('Pump1B', 'Junction1', D=0.15, L=20)

    # Connect to series pump
    G.add_edge('Junction1', 'Pump2', D=0.15, L=50)
    G.add_edge('Pump2', 'Junction2', D=0.15, L=30)
    G.add_edge('Junction2', 'Destination', D=0.15, L=100)

    return G

def analyze_pump_configuration(G):
    """Analyze pump system configuration"""
    pumps = [n for n, d in G.nodes(data=True) if d.get('type') == 'pump']

    print(f"Total pumps in system: {len(pumps)}")
    print(f"Pump nodes: {pumps}")

    # Find pumps in series (on same path)
    # Find pumps in parallel (different paths between same nodes)

    for pump in pumps:
        in_degree = G.in_degree(pump)
        out_degree = G.out_degree(pump)
        print(f"\n{pump}:")
        print(f"  Rated head: {G.nodes[pump].get('head', 'N/A')} m")
        print(f"  Rated flow: {G.nodes[pump].get('flow_rated', 'N/A')} m³/s")
        print(f"  Connections: {in_degree} in, {out_degree} out")
```

## Best Practices

1. **Always validate network connectivity** before performing flow analysis
2. **Use consistent units** for all edge attributes (SI recommended)
3. **Include elevation data** for gravity-driven flow analysis
4. **Model check valves** using directed edges
5. **Consider using MultiDiGraph** for parallel pipes between same junctions
6. **Store results** as node/edge attributes for visualization
7. **Verify conservation of mass** at junctions

## Common Patterns

```python
# Pattern 1: Add multiple parallel pipes
G.add_edge('A', 'B', key='pipe1', D=0.1, L=100)
G.add_edge('A', 'B', key='pipe2', D=0.15, L=100)  # Use MultiDiGraph

# Pattern 2: Iterate over all pipes
for u, v, data in G.edges(data=True):
    process_pipe(u, v, data)

# Pattern 3: Find all junctions (nodes with degree > 2)
junctions = [n for n in G.nodes() if G.degree(n) > 2]

# Pattern 4: Export to GraphML for other tools
nx.write_graphml(G, 'network.graphml')

# Pattern 5: Visualize network layout
pos = nx.spring_layout(G)
# Or use hierarchical layout for better piping visualization
pos = nx.nx_agraph.graphviz_layout(G, prog='dot')
```

## See Also

- `reference.md` - Graph algorithms for hydraulics
- `examples.py` - Complete working examples
- NetworkX documentation for advanced algorithms
