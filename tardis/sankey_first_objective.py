import pandas as pd
from tardis import run_tardis
from tardis.io.configuration.config_reader import Configuration

print("Step 1: Setting up TARDIS Configuration...")

try:
    # 1. Loading the blueprint (.yml file) into Python
    config = Configuration.from_yaml("tardis_example.yml")
    
    # 2. Explicitly setting tracking to 'True' (No YAML edit needed!)
    config.montecarlo.tracking.track_rpacket = True
    print("Tracking explicitly enabled in code!")
    
    print("\nStep 2: Running Simulation")
    # 3. Running the simulation with the new modified config
    sim = run_tardis(
        config, 
        log_level="INFO",
        show_progress_bars=False
    )
    print("\nSimulation Complete!")

    print("\nStep 3: Extracting Tracking DataFrame...")
    
    # 4. Extracting the DataFrame (added checks for different TARDIS versions)
    if hasattr(sim, 'transport_state') and getattr(sim.transport_state, 'tracker_full_df', None) is not None:
        df = sim.transport_state.tracker_full_df
    elif hasattr(sim, 'transport') and getattr(sim.transport.transport_state, 'tracker_full_df', None) is not None:
        df = sim.transport.transport_state.tracker_full_df
    elif hasattr(sim.transport, 'tracker_full_df') and getattr(sim.transport, 'tracker_full_df', None) is not None:
        df = sim.transport.tracker_full_df
    else:
        raise ValueError("DataFrame is still None. Tracking failed to initialize.")
    
    print("\n--- DataFrame Extracted Successfully! ---")
    print("Columns in DataFrame:")
    print(df.columns.tolist())
    
    print("\nFirst 5 rows:")
    # Formatting terminal output so it looks clean
    pd.set_option('display.max_columns', None)
    print(df.head())

    # ---------------------------------------------------------
    # GSoC FIRST OBJECTIVE: Filter and Print Packet History
    # ---------------------------------------------------------
    print("\n" + "="*50)
    print("🚀 EXECUTING GSoC FIRST OBJECTIVE 🚀")
    print("="*50)

    # 1. Filter out 'BOUNDARY' interactions
    df_filtered = df[df['interaction_type'] != 'BOUNDARY']

    # 2. Find the most interesting packet (the one with the most interactions)
    packet_counts = df_filtered.groupby(level='packet_id').size()
    
    if packet_counts.empty:
        print("Oops, no non-boundary interactions found in this run!")
    else:
        best_packet_id = packet_counts.idxmax()
        num_interactions = packet_counts.max()
        
        print(f"\n🎯 Selected Packet ID: {best_packet_id}")
        print(f"🔄 Total non-boundary interactions: {num_interactions}")
        
        # 3. Extract the history for this specific packet
        packet_history = df_filtered.xs(best_packet_id, level='packet_id')
        
        # 4. Print out the list of interactions cleanly
        print("\n📝 Interaction History (Filtered):")
        
        # Selecting important columns to show a clean summary
        cols_to_show = ['interaction_type', 'status', 'radius', 'line_absorb_id', 'line_emit_id']
        clean_history = packet_history[cols_to_show]
        
        print(clean_history.to_string())
        
        # 5. Save to a file to upload to GitHub PR as proof
        clean_history.to_csv("gsoc_packet_history.csv")
        print("\nProof saved to 'gsoc_packet_history.csv'!")
    # ---------------------------------------------------------

except Exception as e:
    print("\nError aagaya:")
    print(e)