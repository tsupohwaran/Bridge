import osmnx as ox
import pandana as pdna
import geopandas as gpd
import pandas as pd
import numpy as np
import os
import warnings
from pathlib import Path

# ================= 1. Global Configuration =================
# Path Configuration

Current_Dir = Path.cwd()
GEO_RAW_DIR = Current_Dir / "data/geo/raw"
GEO_TEMP_DIR = Current_Dir / "data/geo/temp"
GEO_PROCESSED_DIR = Current_Dir / "data/geo/processed"

WORK_DIR = GEO_RAW_DIR
RAW_GRAPH_NAME = "osm_qingdao.graphml"  # Clean road network (no speed, no modification)

# Selected road segment ID file (exported from QGIS as CSV)
BRIDGE_ID_FILE = GEO_RAW_DIR / "bridge_selected.csv"

# Input data
GOV_PATH = GEO_RAW_DIR / "govern_qingdao.xlsx"

# Output files
FIRM_SCENARIOS = [
    {
        "name": "nosec",
        "firm_path": GEO_TEMP_DIR / "firm_qingdao.xlsx",
        "with_bridge": GEO_PROCESSED_DIR / "commute_time_qingdao_python.csv",
        "no_bridge": GEO_PROCESSED_DIR / "commute_time_no_bridge_qingdao_python.csv",
    },
    {
        "name": "model",
        "firm_path": GEO_TEMP_DIR / "firm_list_qingdao_model.xlsx",
        "with_bridge": GEO_PROCESSED_DIR / "commute_time_qingdao_model_python.csv",
        "no_bridge": GEO_PROCESSED_DIR / "commute_time_no_bridge_qingdao_model_python.csv",
    },
    {
        "name": "07_20",
        "firm_path": GEO_RAW_DIR / "firm_qingdao_07_20.xlsx",
        "with_bridge": GEO_PROCESSED_DIR / "commute_time_qingdao_07_20_python.csv",
        "no_bridge": GEO_PROCESSED_DIR / "commute_time_no_bridge_qingdao_07_20_python.csv",
    },
]

# Location
PLACE_NAME = "Qingdao, China"

warnings.filterwarnings('ignore')

# ================= 2. Function Modules =================

def step1_download_raw_network():
    """
    Step 1: Download and save the original road network only.
    No speed assignment, no time calculation, keep it pristine.
    """
    file_path = WORK_DIR / RAW_GRAPH_NAME
    
    if file_path.exists():
        print(f"1. Local road network file detected: {file_path}")
        print("   Skipping download step.")
        return

    print(f"1. Local file not found, downloading road network for '{PLACE_NAME}'...")
    
    # Download drive (vehicle) road network
    G = ox.graph_from_place(PLACE_NAME, network_type='drive')
    
    # Do not add edge speeds at this point, save directly
    print(f"   Saving original road network to: {file_path}")
    ox.save_graphml(G, file_path)
    print("   ✅ Download complete.")

def apply_qingdao_speeds(G):
    """
    Assign speeds to OSM road network based on "Qingdao Central Urban Area Road Network Planning"
    Handles list-type highway tags and non-vehicle roads
    """
    print("   >>> Applying Qingdao Central Urban Area Road Network Planning speed standards...")
    
    # 1. Define base speed mapping (unit: km/h)
    # Arranged from high to low according to OSM tag hierarchy
    speed_lookup = {
        # --- Tier 1: Expressways (58.80) ---
        'motorway': 58.80,      'motorway_link': 58.80,
        'trunk': 58.80,         'trunk_link': 58.80,    # Urban expressways belong here
        
        # --- Tier 2: Arterial roads (34.66) ---
        'primary': 34.66,       'primary_link': 34.66,
        'busway': 34.66,        # Bus-only lanes are usually high-grade
        
        # --- Tier 3: Secondary roads (25.75) ---
        'secondary': 25.75,     'secondary_link': 25.75,
        
        # --- Tier 4: Branch roads (15.82) ---
        'tertiary': 15.82,      'tertiary_link': 15.82,
        'unclassified': 15.82,  # Unclassified roads are usually treated as branch roads
        'highway': 15.82,       # Abnormal generic tag, treat as branch road
        
        # --- Tier 5: Terminal roads/Residential areas (10.00) ---
        'residential': 10.00, 
        'living_street': 10.00,
        
        # --- Special: Non-vehicle (0.001) ---
        'ladder': 0.001         # Stairs, strictly prohibited for vehicles
    }
    
    # Define hierarchy priority (for handling lists, lower number = higher priority)
    hierarchy = {
        'motorway': 1, 'trunk': 2, 
        'primary': 3, 'busway': 3,
        'secondary': 4, 
        'tertiary': 5, 'unclassified': 5,
        'residential': 6, 'living_street': 6
    }

    # 2. Traverse and modify
    count_modified = 0
    for _, _, data in G.edges(data=True):
        highway = data.get('highway')
        
        target_type = None
        
        # --- Case A: Tag is a list (e.g., ['primary', 'trunk']) ---
        if isinstance(highway, list):
            # Find the highest-grade tag in the list (lowest priority number)
            # e.g., ['primary'(3), 'trunk'(2)] -> select 'trunk'
            best_type = None
            min_rank = 999
            
            for h_type in highway:
                # Remove possible link suffix for hierarchy judgment (trunk_link -> trunk)
                base_type = h_type.replace('_link', '')
                rank = hierarchy.get(base_type, 10) # Default to 10 if not found
                
                if rank < min_rank:
                    min_rank = rank
                    best_type = h_type # Keep original tag (including link)
            
            target_type = best_type
            
        # --- Case B: Tag is a string ---
        else:
            target_type = highway
            
        # --- Look up table and assign speed ---
        # If not found in dictionary, default to 10 (residential area speed)
        new_speed = speed_lookup.get(target_type, 10.0)
        
        # Write to maxspeed (used by OSMnx) and speed_kph (custom view)
        data['speed_kph'] = new_speed
        
        # If it's a ladder (stairs), we can even mark it as impassable, but giving very low speed is safest
        if target_type == 'ladder':
            data['speed_kph'] = 0.001
            
        count_modified += 1
        
    print(f"   >>> Updated speed attributes for {count_modified} road segments.")
    return G

def load_target_ids():
    """
    Read QGIS-exported CSV and get all target segment OSMIDs
    """
    df = pd.read_csv(str(BRIDGE_ID_FILE))
    
    # Extract all IDs into a set
    target_ids = set()
    
    for val in df['osmid']:
        # In QGIS-exported CSV, if it's a list, it becomes a string "['123', '456']"
        if isinstance(val, str) and ('[' in val or ',' in val):
            try:
                # Try to parse list string
                # Handle possible non-standard formats, e.g., "123, 456" or "['123', '456']"
                cleaned = val.replace('[', '').replace(']', '').replace("'", "")
                parts = [int(x.strip()) for x in cleaned.split(',') if x.strip().isdigit()]
                target_ids.update(parts)
            except:
                pass
        elif isinstance(val, (int, float)):
            target_ids.add(int(val))
            
    print(f"   Loaded {len(target_ids)} target OSMIDs for blocking.")
    return target_ids

def prepare_scenarios_by_id(edges, target_ids):
    """
    Use ID exact matching to set up experimental scenarios
    (Final version: Fingerprint-level exact matching - Ensure Exact Match)
    """
    print("3. Identifying segments to block based on ID (fingerprint exact matching mode)...")
    
    # 1. Re-read CSV, this time we treat it as "string fingerprints"
    # Even if it contains numbers, we force convert to string to ensure format consistency with edges
    if BRIDGE_ID_FILE.exists():
        df_csv = pd.read_csv(str(BRIDGE_ID_FILE))
        # Convert all osmid column in CSV to strings, remove possible spaces
        # e.g.: "[123, 456]" or "123"
        target_signatures = set(df_csv['osmid'].astype(str).str.strip().values)
        print(f"   >>> Loaded {len(target_signatures)} target fingerprints.")
    else:
        raise FileNotFoundError(f"Cannot find {BRIDGE_ID_FILE}")

    # 2. Convert osmid in road network data to string fingerprints as well
    # Python list [123, 456] converted to string is also "[123, 456]"
    # This allows 1:1 comparison with CSV
    edges_signatures = edges['osmid'].astype(str).str.strip()
    
    # 3. Perform vectorized comparison (very fast and precise)
    # isin returns True only when strings are exactly equal
    mask_bridge = edges_signatures.isin(target_signatures)
    
    bridge_count = mask_bridge.sum()
    print(f"   >>> Precisely identified {bridge_count} road segments as bridges/tunnels.")

    # ================= Export verification file =================
    if bridge_count > 0:
        verify_file = Current_Dir / "output/logs/check_blocked_segments_exact.csv"
        print(f"   >>> Exporting these {bridge_count} segments to: {verify_file}")
        
        blocked_df = edges[mask_bridge].copy()
        
        # Organize columns
        cols = ['name', 'highway', 'osmid', 'length', 'speed_kph', 'u', 'v']
        cols = [c for c in cols if c in blocked_df.columns] 
        other_cols = [c for c in blocked_df.columns if c not in cols and c != 'geometry']
        
        blocked_df[cols + other_cols].to_csv(str(verify_file), index=False)
        print("   ✅ Verification file exported. Please check: the count should now exactly match what you selected in QGIS.")
    else:
        print("   ⚠️ Warning: No segments identified!")
        print("   Please check if the format in CSV is slightly different from Python (e.g., spaces, quotes, etc.).")
        # Print a sample for debugging
        print(f"   CSV sample: {list(target_signatures)[0]}")
        print(f"   Edge sample: {edges_signatures.iloc[0]}")
    # =======================================================

    # Scenario A: Real world (with bridges)
    edges['time_scenario_A'] = edges['travel_time'].copy()
    
    # Scenario B: No bridges - implemented by removing edges, not by giving large weights
    # (pandana's contraction hierarchy has issues with extreme weight differences)
    
    # ================= Diagnostic check =================
    print("\n   >>> 🔍 Diagnostic check:")
    print(f"   Number of blocked segments: {mask_bridge.sum()}")
    
    if mask_bridge.sum() > 0:
        blocked_times_A = edges.loc[mask_bridge, 'time_scenario_A']
        
        print(f"   Scenario A (with bridges) marked segment travel times: min={blocked_times_A.min():.2f}, max={blocked_times_A.max():.2f}, mean={blocked_times_A.mean():.2f}")
        print("   Scenario B (no bridges) will completely remove these segments")
            
        # Check global statistics
        print(f"\n   Global comparison:")
        print(f"   Scenario A edge count: {len(edges)}")
        print(f"   Scenario B edge count: {len(edges) - mask_bridge.sum()}")
    # ===========================================
    
    return edges, mask_bridge

def run_pandana_calc(nodes, edges, time_col, firms, govs, outfile, remove_mask=None):
    """
    Calculate shortest path matrix
    
    Parameters:
    -----------
    remove_mask : pd.Series or None
        Boolean mask marking edges to remove from network (for "no bridge" scenario)
    """
    print(f"\n>>> Calculating matrix: [{time_col}]")
    
    # ===== Key fix: Correctly handle node ID to index mapping =====
    # pandana requires u, v to be node index positions in the array (0, 1, 2, ...)
    # But OSMnx gives u, v as actual OSM node IDs (very large numbers)
    
    # 1. Ensure nodes have x, y columns
    nodes_clean = nodes.copy()
    if 'x' not in nodes_clean.columns:
        nodes_clean['x'] = nodes_clean.geometry.x
        nodes_clean['y'] = nodes_clean.geometry.y
    
    # 2. Create OSM ID -> index position mapping
    # nodes index is osmid
    node_ids = nodes_clean.index.values
    node_id_to_idx = {nid: idx for idx, nid in enumerate(node_ids)}
    
    # 3. Prepare edge data
    # If remove_mask exists, remove these edges (instead of giving large weights)
    if remove_mask is not None:
        edges_to_use = edges[~remove_mask].copy()
        print(f"   Removed {remove_mask.sum()} edges (bridges/tunnels)")
    else:
        edges_to_use = edges.copy()
    
    edges_clean = pd.DataFrame({
        'u': edges_to_use['u'].values,
        'v': edges_to_use['v'].values,
        'weight': edges_to_use[time_col].values
    })
    
    edges_clean['u_idx'] = edges_clean['u'].map(node_id_to_idx)
    edges_clean['v_idx'] = edges_clean['v'].map(node_id_to_idx)
    
    # Check for unmappable nodes (should not happen theoretically)
    missing = edges_clean['u_idx'].isna().sum() + edges_clean['v_idx'].isna().sum()
    if missing > 0:
        print(f"   ⚠️ Warning: {missing} edge node IDs could not be mapped!")
        edges_clean = edges_clean.dropna(subset=['u_idx', 'v_idx'])
    
    print(f"   Node count: {len(nodes_clean)}, Edge count: {len(edges_clean)}")
    print(f"   {time_col} stats: min={edges_clean['weight'].min():.2f}, max={edges_clean['weight'].max():.2f}, mean={edges_clean['weight'].mean():.2f}")
    
    # 4. Create pandana network
    net = pdna.Network(
        nodes_clean['x'].values, 
        nodes_clean['y'].values,
        edges_clean['u_idx'].astype(int).values, 
        edges_clean['v_idx'].astype(int).values,
        edges_clean[['weight']],
        twoway=True 
    )
    
    # Coordinate system alignment
    firms_proj = firms.to_crs("EPSG:4326") if firms.crs.to_string() != "EPSG:4326" else firms
    govs_proj = govs.to_crs("EPSG:4326") if govs.crs.to_string() != "EPSG:4326" else govs
    
    firm_nodes = net.get_node_ids(firms_proj.geometry.x.values, firms_proj.geometry.y.values)
    gov_nodes = net.get_node_ids(govs_proj.geometry.x.values, govs_proj.geometry.y.values)
    
    sources = np.tile(firm_nodes, len(gov_nodes))
    targets = np.repeat(gov_nodes, len(firm_nodes))
    
    print("   Starting high-speed calculation...")
    times = net.shortest_path_lengths(sources, targets)
    
    df = pd.DataFrame({
        'firm_id': np.tile(firms_proj['id'].values, len(gov_nodes)),
        'town_id': np.repeat(govs_proj['town'].values, len(firm_nodes)),
        'travel_time_sec': times
    })
    
    df.loc[df['travel_time_sec'] > 1e6, 'travel_time_sec'] = np.nan
    df['travel_time_min'] = (df['travel_time_sec'] / 60).round(2)
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(outfile), index=False)
    print(f"   ✅ Complete: {outfile}")

# ================= 3. Main Process =================

def main():
    # --- Step 1: Download (skip if already exists) ---
    step1_download_raw_network()
    
    # --- Step 2: Load and process ---
    print("\n2. Loading road network and applying Qingdao speed standards...")
    file_path = WORK_DIR / RAW_GRAPH_NAME
    G = ox.load_graphml(file_path)
    
    # Apply Qingdao speeds
    G = apply_qingdao_speeds(G)
    
    # Calculate time weights (Time = Length / Speed)
    G = ox.add_edge_travel_times(G)
    
    # Convert to dataframes
    nodes, edges = ox.graph_to_gdfs(G, nodes=True, edges=True)
    edges = edges.reset_index()
    
    # --- Step 3: Set up scenarios (ID exact matching) ---
    # Load IDs exported from QGIS
    target_ids = load_target_ids()
    
    # Mark bridges, return edges and mask
    edges, bridge_mask = prepare_scenarios_by_id(edges, target_ids)
    
    # --- Step 4: Run calculations ---
    df_govs = pd.read_excel(str(GOV_PATH))
    gdf_govs = gpd.GeoDataFrame(
        df_govs,
        geometry=gpd.points_from_xy(df_govs['longitude'], df_govs['latitude']),
        crs="EPSG:4326"
    )

    for scenario in FIRM_SCENARIOS:
        if not scenario["firm_path"].exists():
            print(f"\n4. Skipping firm sample '{scenario['name']}' because {scenario['firm_path']} does not exist.")
            continue

        print(f"\n4. Running firm sample: {scenario['name']}")
        df_firms = pd.read_excel(str(scenario["firm_path"]))
        gdf_firms = gpd.GeoDataFrame(
            df_firms,
            geometry=gpd.points_from_xy(df_firms['longitude'], df_firms['latitude']),
            crs="EPSG:4326"
        )

        # Scenario A: With bridges (use complete road network)
        run_pandana_calc(nodes, edges, 'time_scenario_A', gdf_firms, gdf_govs, scenario["with_bridge"], remove_mask=None)

        # Scenario B: No bridges (remove bridge/tunnel road segments)
        run_pandana_calc(nodes, edges, 'time_scenario_A', gdf_firms, gdf_govs, scenario["no_bridge"], remove_mask=bridge_mask)

if __name__ == "__main__":
    main()
