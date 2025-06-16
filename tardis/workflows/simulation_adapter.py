class SimulationAdapter:
    def __init__(self, workflow):
        self._workflow = workflow
    
    def __getattr__(self, name):
        return getattr(self._workflow, name)
    
    @property
    def plasma(self):
        return self._workflow.plasma_solver
    
    @property 
    def transport(self):
        return self._workflow.transport_solver
        
    @property
    def iterations_executed(self):
        return self._workflow.completed_iterations
        
    @property
    def iterations(self):
        return self._workflow.total_iterations
        
    @property
    def no_of_packets(self):
        return self._workflow.real_packet_count
        
    @property
    def last_no_of_packets(self):
        return self._workflow.final_iteration_packet_count
        
    @property
    def no_of_virtual_packets(self):
        return self._workflow.virtual_packet_count


