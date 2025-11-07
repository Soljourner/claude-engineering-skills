"""
NetworkX Flow Networks Examples
Hydraulic piping systems and multi-pump configurations
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# Example 1: Simple Piping Network Model
# ============================================================================

def example_1_simple_piping_network():
    """Create and analyze a simple piping network"""
    print("=" * 70)
    print("EXAMPLE 1: Simple Piping Network")
    print("=" * 70)

    # Create directed graph
    G = nx.DiGraph()

    # Add nodes with properties
    nodes_data = {
        'Reservoir': {'type': 'source', 'elevation': 50.0, 'pressure': 101325},
        'Junction1': {'type': 'junction', 'elevation': 10.0},
        'Junction2': {'type': 'junction', 'elevation': 5.0},
        'Junction3': {'type': 'junction', 'elevation': 0.0},
        'Outlet1': {'type': 'sink', 'elevation': 0.0, 'demand': 0.02},
        'Outlet2': {'type': 'sink', 'elevation': 0.0, 'demand': 0.03}
    }

    for node, attrs in nodes_data.items():
        G.add_node(node, **attrs)

    # Add pipes with hydraulic properties
    # Format: (from, to, diameter_m, length_m, roughness_m)
    pipes_data = [
        ('Reservoir', 'Junction1', 0.20, 100, 0.00015),
        ('Junction1', 'Junction2', 0.15, 80, 0.00015),
        ('Junction1', 'Junction3', 0.10, 60, 0.00015),
        ('Junction2', 'Outlet1', 0.10, 50, 0.00015),
        ('Junction3', 'Outlet2', 0.10, 70, 0.00015),
    ]

    for u, v, D, L, k in pipes_data:
        # Calculate pipe resistance coefficient (K = f*L/D, assume f=0.02)
        f = 0.02
        K = f * L / D
        # Calculate cross-sectional area
        A = np.pi * (D/2)**2

        G.add_edge(u, v, D=D, L=L, k=k, K=K, A=A)

    # Network statistics
    print(f"\nNetwork Statistics:")
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Pipes: {G.number_of_edges()}")
    print(f"  Average node degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")

    # Check connectivity
    is_connected = nx.is_weakly_connected(G)
    print(f"  Network is connected: {is_connected}")

    # Find all paths from reservoir to outlets
    print(f"\nFlow Paths:")
    outlets = ['Outlet1', 'Outlet2']
    for outlet in outlets:
        paths = list(nx.all_simple_paths(G, 'Reservoir', outlet))
        print(f"\n  To {outlet} ({len(paths)} paths):")
        for i, path in enumerate(paths, 1):
            total_length = sum(G[u][v]['L'] for u, v in zip(path[:-1], path[1:]))
            min_diameter = min(G[u][v]['D'] for u, v in zip(path[:-1], path[1:]))
            print(f"    Path {i}: {' -> '.join(path)}")
            print(f"            Length: {total_length} m, Min diameter: {min_diameter} m")

    return G


# ============================================================================
# Example 2: Multi-Pump System Graph
# ============================================================================

def example_2_multipump_system():
    """Model a multi-pump system with parallel and series configurations"""
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Multi-Pump System")
    print("=" * 70)

    G = nx.DiGraph()

    # System components
    components = {
        'Tank_A': {'type': 'tank', 'elevation': 0, 'volume': 50},
        'Suction_Header': {'type': 'junction', 'elevation': 0},
        'Pump_1A': {'type': 'pump', 'head': 45, 'flow_rated': 0.04, 'power': 15000},
        'Pump_1B': {'type': 'pump', 'head': 45, 'flow_rated': 0.04, 'power': 15000},
        'Discharge_Header': {'type': 'junction', 'elevation': 0},
        'Booster_Pump': {'type': 'pump', 'head': 30, 'flow_rated': 0.08, 'power': 20000},
        'Elevated_Tank': {'type': 'tank', 'elevation': 60, 'volume': 100}
    }

    for node, attrs in components.items():
        G.add_node(node, **attrs)

    # Connections
    connections = [
        ('Tank_A', 'Suction_Header', 0.20, 15),
        ('Suction_Header', 'Pump_1A', 0.15, 5),
        ('Suction_Header', 'Pump_1B', 0.15, 5),
        ('Pump_1A', 'Discharge_Header', 0.15, 10),
        ('Pump_1B', 'Discharge_Header', 0.15, 10),
        ('Discharge_Header', 'Booster_Pump', 0.15, 50),
        ('Booster_Pump', 'Elevated_Tank', 0.15, 120),
    ]

    for u, v, D, L in connections:
        A = np.pi * (D/2)**2
        G.add_edge(u, v, D=D, L=L, A=A)

    # Analyze pump configuration
    print(f"\nPump Configuration Analysis:")
    pumps = [n for n, d in G.nodes(data=True) if d.get('type') == 'pump']
    print(f"  Total pumps: {len(pumps)}")

    for pump in pumps:
        data = G.nodes[pump]
        print(f"\n  {pump}:")
        print(f"    Head: {data['head']} m")
        print(f"    Rated flow: {data['flow_rated']} m³/s ({data['flow_rated']*3600:.1f} m³/h)")
        print(f"    Power: {data['power']/1000:.1f} kW")

        # Check configuration
        predecessors = list(G.predecessors(pump))
        successors = list(G.successors(pump))
        print(f"    Input from: {predecessors}")
        print(f"    Output to: {successors}")

    # Identify parallel pumps (share common predecessor and successor)
    print(f"\n  Parallel Pump Groups:")
    parallel_pumps = [['Pump_1A', 'Pump_1B']]
    for group in parallel_pumps:
        total_capacity = sum(G.nodes[p]['flow_rated'] for p in group)
        total_power = sum(G.nodes[p]['power'] for p in group)
        print(f"    {' & '.join(group)}")
        print(f"      Combined capacity: {total_capacity} m³/s ({total_capacity*3600:.1f} m³/h)")
        print(f"      Combined power: {total_power/1000:.1f} kW")

    # Calculate total system head
    pump_heads = [G.nodes[p]['head'] for p in pumps]
    print(f"\n  System Characteristics:")
    print(f"    Total head (series pumps): {sum([G.nodes['Booster_Pump']['head']])} m")
    print(f"    Parallel stage head: {G.nodes['Pump_1A']['head']} m")
    print(f"    Total elevation gain: {G.nodes['Elevated_Tank']['elevation']} m")

    return G


# ============================================================================
# Example 3: Flow Distribution Calculation
# ============================================================================

def example_3_flow_distribution():
    """Calculate flow distribution using max flow algorithm"""
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Flow Distribution Analysis")
    print("=" * 70)

    # Create a branched network
    G = nx.DiGraph()

    nodes = {
        'Source': {'type': 'source', 'elevation': 30},
        'J1': {'type': 'junction', 'elevation': 20},
        'J2': {'type': 'junction', 'elevation': 15},
        'J3': {'type': 'junction', 'elevation': 15},
        'J4': {'type': 'junction', 'elevation': 10},
        'Sink1': {'type': 'sink', 'elevation': 0, 'demand': 0.015},
        'Sink2': {'type': 'sink', 'elevation': 0, 'demand': 0.025},
    }

    for node, attrs in nodes.items():
        G.add_node(node, **attrs)

    # Pipes with different diameters (creates unequal flow distribution)
    pipes = [
        ('Source', 'J1', 0.20, 50),
        ('J1', 'J2', 0.15, 80),
        ('J1', 'J3', 0.10, 60),  # Smaller diameter path
        ('J2', 'J4', 0.15, 70),
        ('J3', 'J4', 0.10, 50),
        ('J4', 'Sink1', 0.10, 40),
        ('J2', 'Sink2', 0.15, 90),
    ]

    for u, v, D, L in pipes:
        A = np.pi * (D/2)**2
        # Capacity based on reasonable velocity limit (2.5 m/s)
        capacity = A * 2.5
        G.add_edge(u, v, D=D, L=L, A=A, capacity=capacity)

    # Perform max flow analysis
    print(f"\nNetwork Configuration:")
    print(f"  Pipes: {G.number_of_edges()}")
    print(f"  Junctions: {G.number_of_nodes()}")

    # Calculate max flow from source to each sink
    sinks = ['Sink1', 'Sink2']

    print(f"\nMax Flow Analysis:")
    for sink in sinks:
        flow_value, flow_dict = nx.maximum_flow(G, 'Source', sink)
        demand = G.nodes[sink]['demand']

        print(f"\n  To {sink}:")
        print(f"    Maximum capacity: {flow_value:.4f} m³/s ({flow_value*3600:.1f} m³/h)")
        print(f"    Required demand: {demand:.4f} m³/s ({demand*3600:.1f} m³/h)")
        print(f"    Adequacy: {'OK' if flow_value >= demand else 'INSUFFICIENT'}")

        print(f"    Flow distribution:")
        for u in flow_dict:
            for v, flow in flow_dict[u].items():
                if flow > 0:
                    edge_data = G[u][v]
                    velocity = flow / edge_data['A']
                    utilization = flow / edge_data['capacity'] * 100
                    print(f"      {u} -> {v}: {flow:.4f} m³/s, "
                          f"v={velocity:.2f} m/s, util={utilization:.1f}%")

    # Find minimum cut (bottleneck)
    print(f"\nBottleneck Analysis:")
    for sink in sinks:
        cut_value, partition = nx.minimum_cut(G, 'Source', sink)
        reachable, non_reachable = partition

        print(f"\n  To {sink}:")
        print(f"    Min cut capacity: {cut_value:.4f} m³/s")
        print(f"    Bottleneck pipes:")
        for u, v, data in G.edges(data=True):
            if u in reachable and v in non_reachable:
                print(f"      {u} -> {v}: D={data['D']*1000:.0f} mm, "
                      f"capacity={data['capacity']:.4f} m³/s")

    return G


# ============================================================================
# Example 4: Network Visualization
# ============================================================================

def example_4_network_visualization():
    """Visualize piping network with flow data"""
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Network Visualization")
    print("=" * 70)

    # Create a simple network
    G = nx.DiGraph()

    # Add nodes
    nodes = ['Tank', 'P1', 'J1', 'J2', 'Out1', 'Out2']
    for i, node in enumerate(nodes):
        if node.startswith('Tank'):
            node_type = 'source'
        elif node.startswith('P'):
            node_type = 'pump'
        elif node.startswith('J'):
            node_type = 'junction'
        else:
            node_type = 'sink'
        G.add_node(node, type=node_type)

    # Add edges
    edges = [
        ('Tank', 'P1', 0.15, 20),
        ('P1', 'J1', 0.15, 50),
        ('J1', 'J2', 0.10, 40),
        ('J1', 'Out1', 0.10, 30),
        ('J2', 'Out2', 0.10, 35),
    ]

    for u, v, D, L in edges:
        A = np.pi * (D/2)**2
        G.add_edge(u, v, D=D, L=L, A=A, capacity=A*2.0)

    # Calculate flow distribution
    flow_value, flow_dict = nx.maximum_flow(G, 'Tank', 'Out1')

    # Add flow to edges
    for u in flow_dict:
        for v, flow in flow_dict[u].items():
            if G.has_edge(u, v):
                G[u][v]['flow'] = flow

    print(f"\nGenerating visualization...")

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Layout - hierarchical left to right
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Plot 1: Network topology with pipe diameters
    ax1.set_title('Piping Network Topology', fontsize=14, fontweight='bold')

    # Node colors by type
    node_colors = []
    for node in G.nodes():
        node_type = G.nodes[node]['type']
        if node_type == 'source':
            node_colors.append('lightblue')
        elif node_type == 'pump':
            node_colors.append('orange')
        elif node_type == 'junction':
            node_colors.append('lightgray')
        else:
            node_colors.append('lightgreen')

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000,
                           ax=ax1, edgecolors='black', linewidths=2)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax1)

    # Edge widths proportional to diameter
    edge_widths = [G[u][v]['D'] * 20 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='gray',
                          arrows=True, arrowsize=20, ax=ax1,
                          connectionstyle='arc3,rad=0.1')

    # Edge labels with diameter
    edge_labels = {(u, v): f"D={d['D']*1000:.0f}mm\nL={d['L']:.0f}m"
                   for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=7, ax=ax1)

    ax1.axis('off')
    ax1.set_aspect('equal')

    # Plot 2: Flow distribution
    ax2.set_title('Flow Distribution', fontsize=14, fontweight='bold')

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=1000,
                          ax=ax2, edgecolors='black', linewidths=2)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax2)

    # Edge colors by flow
    edges_with_flow = [(u, v) for u, v in G.edges() if G[u][v].get('flow', 0) > 0]
    flows = [G[u][v].get('flow', 0) for u, v in edges_with_flow]

    if flows:
        max_flow = max(flows)
        edge_colors = [G[u][v].get('flow', 0) / max_flow for u, v in G.edges()]
        edge_widths = [5 + G[u][v].get('flow', 0) / max_flow * 10 for u, v in G.edges()]

        nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color=edge_colors,
                              edge_cmap=plt.cm.Reds, arrows=True, arrowsize=20,
                              ax=ax2, connectionstyle='arc3,rad=0.1')

        # Edge labels with flow
        edge_labels = {(u, v): f"{d.get('flow', 0)*3600:.1f} m³/h"
                       for u, v, d in G.edges(data=True) if d.get('flow', 0) > 0}
        nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, ax=ax2)

    ax2.axis('off')
    ax2.set_aspect('equal')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', label='Source/Tank'),
        Patch(facecolor='orange', edgecolor='black', label='Pump'),
        Patch(facecolor='lightgray', edgecolor='black', label='Junction'),
        Patch(facecolor='lightgreen', edgecolor='black', label='Sink/Outlet')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, frameon=True)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)

    # Save figure
    output_file = '/tmp/networkx_flow_network.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  Visualization saved to: {output_file}")

    # Display if in interactive mode
    try:
        plt.show()
    except:
        print(f"  (Non-interactive environment - plot saved to file)")

    return G


# ============================================================================
# Main execution
# ============================================================================

if __name__ == '__main__':
    print("\n")
    print("=" * 70)
    print(" NetworkX Flow Networks - Hydraulic Systems Examples")
    print("=" * 70)

    # Run all examples
    G1 = example_1_simple_piping_network()
    G2 = example_2_multipump_system()
    G3 = example_3_flow_distribution()
    G4 = example_4_network_visualization()

    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
    print("\nKey Takeaways:")
    print("  - NetworkX provides powerful graph algorithms for flow networks")
    print("  - Piping systems map naturally to directed graphs")
    print("  - Max flow algorithms identify system capacity and bottlenecks")
    print("  - Path analysis helps optimize network design")
    print("  - Visualization aids in understanding complex networks")
    print("=" * 70 + "\n")
