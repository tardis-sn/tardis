.. .. pyodide::
..     from pyodide.http import pyfetch
..     import panel as pn
..     from shell_info import ShellInfoWidget, HDFShellInfo

..     pn.extension('tabulator')

..     async def fetch_and_load():
..         url = "https://archive.org/download/sim_20250329/sim.hdf"
..         response = await pyfetch(url)
        
..         with open("sim.hdf", "wb") as f:
..             f.write(await response.bytes())

..         shell_info_data = HDFShellInfo("sim.hdf")
..         shell_info_widget = ShellInfoWidget(shell_info_data)
..         return shell_info_widget.display()

..     app = await fetch_and_load()
..     app