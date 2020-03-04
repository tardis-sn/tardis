'''
This is a script for easily configuring a local dask cluster.

Run this before running the "run_parallel_sims.py" file in another terminal.
'''
import subprocess

def start_local_cluster(hostname,port,n_processes,n_threads):
    '''
    Start a local dask cluster.
    -- hostname TEXT         - hostname for scheduler
    -- port TEXT     	     - port number for scheduler
    -- n_processes INTEGER   - number of worker nodes
    -- n_threads INTEGER     - number of threads per worker

    * Example
        >> start_local_cluster(hostname='192.xxx.x.x',port='5500',n_processes=20,n_threads=3)
        >> # this will launch a dask-scheduler along with 20 dask-workers registered to it each with a ThreadPoolExecutor of size 3.
    '''
    #since the process of launching a scheduler and the worker nodes is parallel, some warnings may show up in the terminal
    subprocess.Popen(f'dask-scheduler --host {hostname} --port {port}',shell=True)
    subprocess.Popen(f'dask-worker {hostname}:{port} --nprocs {str(n_processes)} --nthreads {str(n_threads)}',shell=True)

def main():
	hostname=input('Scheduler Hostname (run "dask-scheduler" in terminal to see details) :> ')
	port=input('Scheduler Port Number :> ')
	n_processes=int(input('Number of worker nodes :> '))
	n_threads=int(input('Number of threads per node :> '))

	start_local_cluster(hostname=hostname,port=port,n_processes=n_processes,n_threads=n_threads)

if __name__=='__main__':
	main()
