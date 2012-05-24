import plasma
import initialize

def run_oned(conn, fname):
    initial_config = initialize.read_simple_config(fname)
    return initial_config