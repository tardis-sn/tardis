import pandas as pd
import matplotlib.pyplot as plt

def run_sankey_logic():
    # 1. Representative Data Structure
    data = {
        'packet_id': [22, 22, 22, 22, 45, 45, 45, 45, 45],
        'interaction_type': [
            'boundary', 'electron_scattering', 'boundary', 'line_interaction',
            'boundary', 'electron_scattering', 'line_interaction', 'line_interaction', 'escape'
        ],
        'r': [1.1e15, 1.2e15, 1.3e15, 1.4e15,
              1.1e15, 1.2e15, 1.3e15, 1.4e15, 1.6e15]
    }
    tracker_full_df = pd.DataFrame(data)

    # 2. Filtering Logic: Remove boundary interactions
    filtered_df = tracker_full_df[tracker_full_df['interaction_type'] != 'boundary'].copy()

    # 3. Identify most active packet automatically
    physical_interactions = filtered_df[filtered_df['interaction_type'] != 'escape']
    best_packet_id = physical_interactions.groupby('packet_id').size().idxmax()

    # 4. Extract history for selected packet
    packet_history = filtered_df[filtered_df['packet_id'] == best_packet_id].reset_index(drop=True)
    
    print(f"=== Packet {best_packet_id} History (Filtered) ===")
    print(packet_history[['packet_id', 'interaction_type', 'r']])

    # 5. Step Plot Visualization
    plt.figure(figsize=(10, 5))
    plt.step(range(len(packet_history)), packet_history['r'], where='post', marker='o', color='purple')

    for i, row in packet_history.iterrows():
        plt.annotate(row['interaction_type'], (i, row['r']), textcoords="offset points", xytext=(0,10), ha='center')

    plt.title(f"Packet {best_packet_id} Interaction Path (First Objective)")
    plt.xlabel("Interaction Step")
    plt.ylabel("Radius (r)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()

if __name__ == "__main__":
    run_sankey_logic()
