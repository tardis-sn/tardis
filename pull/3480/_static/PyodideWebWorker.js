importScripts("https://cdn.jsdelivr.net/pyodide/v0.28.2/full/pyodide.js");

const QUEUE = [];

const REQUIRES = {}

function sendPatch(patch, buffers, cell_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers,
    id: cell_id
  })
}

function sendStdout(cell_id, stdout) {
  self.postMessage({
    type: 'stdout',
    content: stdout,
    id: cell_id
  })
}
function sendStderr(cell_id, stderr) {
  self.postMessage({
    type: 'stderr',
    content: stderr,
    id: cell_id
  })
}

async function loadApplication(cell_id, path) {
  console.log("Loading pyodide!");
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  self.pyodide.globals.set("sendStdout", sendStdout);
  self.pyodide.globals.set("sendStderr", sendStderr);
  console.log("Loaded!");
  await self.pyodide.loadPackage("micropip");
  const packages = ['panel', 'pandas', 'bokeh', 'pyodide-http', 'cloudpickle', 'astropy', 'numpy', 'scipy', 'matplotlib', 'plotly'];
  if (path != null) {
    for (const key of Object.keys(REQUIRES)) {
      if (path.replace('.html', '').endsWith(key.replace('.md', ''))) {
        for (const req of REQUIRES[key]) {
          packages.push(req)
        }
      }
    }
  }

  await self.pyodide.runPythonAsync("\"\"\"\nSetup code to run when any pyodide cell is executed.\nThis is responsible for:\n    1. Creating lightweight mock modules for certain tardis modules so that\n        pyodide does not need to install the whole dependency tree of tardis.\n    2. Writing the files needed by pyodide into the virtual filesystem.\n\nNote:\n    This file is executed in a Pyodide web worker using pyodide.runPythonAsync().\n    During execution, modules such as pyodide.http and js are available.\n    The code runs asynchronously.\n\"\"\"\n\nimport sysconfig\nimport pyodide.http\nimport js\nfrom pathlib import Path\n\n\nlib = Path(sysconfig.get_path(\'purelib\'))  # Path to site-packages in pyodide\'s virtual filesystem\n\n# Temporarily, this points to my fork\n# Can be changed to the main repo later\nGITHUB_URL = \'https://raw.githubusercontent.com/SS-9098/tardis/refs/heads/First-Objective-Migrate-to-Panel/\'\n# Derive _static/ base URL from the WebWorker\'s own script URL.\n# The worker is always loaded from _static/PyodideWebWorker.js.\nWORKER_URL = js.self.location.href\nSTATIC_URL = WORKER_URL[: WORKER_URL.rfind(\'/\') + 1]  # trailing slash\n\ndef write_package(rel_path, content):\n    full = lib / rel_path\n    full.parent.mkdir(parents=True, exist_ok=True)\n    full.write_text(content)\n\n\ndef write_data(rel_path, content):\n    p = Path(rel_path)\n    p.parent.mkdir(parents=True, exist_ok=True)\n    if p.suffix == \'.pkl\':\n        p.write_bytes(content)\n    else:\n        p.write_text(content)\n\n\n# Package __init__.py files\nfor _pkg in [\n    \'tardis/__init__.py\',\n    \'tardis/util/__init__.py\',\n    \'tardis/configuration/__init__.py\',\n    \'tardis/transport/__init__.py\',\n    \'tardis/transport/montecarlo/__init__.py\',\n    \'tardis/transport/montecarlo/packets/__init__.py\',\n    \'tardis/visualization/__init__.py\',\n    \'tardis/visualization/widgets/__init__.py\',\n]:\n    write_package(_pkg, \'\')\n\nwrite_package(\'tardis/visualization/__init__.py\', \"\"\"\\\nfrom tardis.visualization.widgets.shell_info import shell_info_from_simulation\nfrom tardis.visualization.widgets.line_info import LineInfoWidget\nimport panel as pn\npn.extension(\"tabulator\", \"plotly\")\n\"\"\")\n\n# stub of tardis.util.base, which is needed by the widgets\nwrite_package(\'tardis/util/base.py\', \"\"\"\\\nimport pandas as pd\nimport re\nfrom collections import OrderedDict\n\nATOMIC_SYMBOLS_DATA = (\n    pd.read_csv(\n        \"/data/atomic_symbols.dat\",\n        sep=r\"\\s+\",\n        names=[\"atomic_number\", \"symbol\"],\n    )\n    .set_index(\"atomic_number\")\n    .squeeze()\n)\n\nATOMIC_NUMBER2SYMBOL = OrderedDict(ATOMIC_SYMBOLS_DATA.to_dict())\nSYMBOL2ATOMIC_NUMBER = OrderedDict((y, x) for x, y in ATOMIC_NUMBER2SYMBOL.items())\n\nNUMERAL_MAP = tuple(\n    zip(\n        (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1),\n        (\"M\", \"CM\", \"D\", \"CD\", \"C\", \"XC\", \"L\", \"XL\", \"X\", \"IX\", \"V\", \"IV\", \"I\"),\n    )\n)\n\ndef int_to_roman(i):\n    result = []\n    for integer, numeral in NUMERAL_MAP:\n        count = i // integer\n        result.append(numeral * count)\n        i -= integer * count\n    return \"\".join(result)\n\ndef roman_to_int(roman_string):\n    NUMERALS_SET = set(list(zip(*NUMERAL_MAP))[1])\n    roman_string = roman_string.upper()\n    if len(set(roman_string.upper()) - NUMERALS_SET) != 0:\n        raise ValueError(f\"{roman_string} does not seem to be a roman numeral\")\n    i = result = 0\n    for integer, numeral in NUMERAL_MAP:\n        while roman_string[i : i + len(numeral)] == numeral:\n            result += integer\n            i += len(numeral)\n    if result < 1:\n        raise ValueError(f\"Can not interpret Roman Numeral {roman_string}\")\n    return result\n\ndef reformat_element_symbol(element_string):\n    return element_string[0].upper() + element_string[1:].lower()\n\ndef atomic_number2element_symbol(atomic_number):\n    return ATOMIC_NUMBER2SYMBOL[atomic_number]\n    \ndef element_symbol2atomic_number(element_string):\n    reformatted_element_string = reformat_element_symbol(element_string)\n    return SYMBOL2ATOMIC_NUMBER[reformatted_element_string]\n\ndef species_tuple_to_string(species_tuple, roman_numerals=True):\n    atomic_number, ion_number = species_tuple\n    element_symbol = ATOMIC_NUMBER2SYMBOL[atomic_number]\n    if roman_numerals:\n        roman_ion_number = int_to_roman(ion_number + 1)\n        return f\"{str(element_symbol)} {roman_ion_number}\"\n    return f\"{element_symbol} {ion_number:d}\"\n    \ndef species_string_to_tuple(species_string):\n    try:\n        element_symbol, ion_number_string = re.match(\n            r\"^(\\w+)\\s*(\\d+)\", species_string\n        ).groups()\n    except AttributeError:\n        try:\n            element_symbol, ion_number_string = species_string.split()\n        except ValueError:\n            print(\n                f\'Species string \"{species_string}\" is not of format <element_symbol><number>\'\n                f\" (e.g. Fe 2, Fe2, ..)\"\n            )\n\n    atomic_number = element_symbol2atomic_number(element_symbol)\n\n    try:\n        ion_number = roman_to_int(ion_number_string)\n    except ValueError:\n        print(\n            \"Ion Number does not contain a Roman Numeral. Checking for integer value\"\n        )\n        try:\n            ion_number = int(ion_number_string)\n        except ValueError:\n            print(\n                f\"Given ion number (\'{ion_number_string}\') could not be parsed\"\n            )\n\n    if ion_number - 1 > atomic_number:\n        raise print(\"Species given does not exist: ion number > atomic number\")\n\n    return atomic_number, ion_number - 1\n\"\"\")\n\nwrite_package(\'tardis/util/environment.py\', \"\"\"\\\nclass Environment:\n    @staticmethod\n    def allows_widget_display():\n        return True\n\n    @classmethod\n    def get_current_environment(cls):\n        return \'pyodide\'\n\"\"\")\n\nwrite_package(\'tardis/transport/montecarlo/packets/radiative_packet.py\', \"\"\"\\\nfrom enum import IntEnum\nclass InteractionType(IntEnum):\n    NO_INTERACTION = -1\n    BOUNDARY = 1\n    LINE = 2\n    ESCATTERING = 4\n    CONTINUUM_PROCESS = 8\n\"\"\")\n\nwrite_package(\'tardis/configuration/sorting_globals.py\', \"\"\"\\\nSORTING_ALGORITHM = \"stable\"\n\"\"\")\n\nwrite_package(\'tardis/constants.py\', \"\"\"\\\nfrom astropy.constants.astropyconst13 import *\n\"\"\")\n\nresp = await pyodide.http.pyfetch(f\'{GITHUB_URL}tardis/visualization/widgets/util.py\')\nwrite_package(\'tardis/visualization/widgets/util.py\', await resp.string())\n\nresp = await pyodide.http.pyfetch(f\'{GITHUB_URL}tardis/visualization/widgets/shell_info.py\')\nwrite_package(\'tardis/visualization/widgets/shell_info.py\', await resp.string())\n\nresp = await pyodide.http.pyfetch(f\'{GITHUB_URL}tardis/analysis.py\')\nwrite_package(\'tardis/analysis.py\', await resp.string())\n\nresp = await pyodide.http.pyfetch(f\'{GITHUB_URL}tardis/visualization/widgets/line_info.py\')\nwrite_package(\'tardis/visualization/widgets/line_info.py\', await resp.string())\n\nresp = await pyodide.http.pyfetch(f\'{GITHUB_URL}tardis/visualization/plot_util.py\')\nwrite_package(\'tardis/visualization/plot_util.py\', await resp.string())\n\nresp = await pyodide.http.pyfetch(STATIC_URL + \'data/atomic_symbols.dat\')\nwrite_data(\'/data/atomic_symbols.dat\', await resp.string())\n\nresp = await pyodide.http.pyfetch(STATIC_URL + \'data/line_info_data.pkl\')\nwrite_data(\'/data/line_info_data.pkl\', await resp.bytes())\n\nresp = await pyodide.http.pyfetch(STATIC_URL + \'data/shell_info_data.pkl\')\nwrite_data(\'/data/shell_info_data.pkl\', await resp.bytes())\n")
  self.pyodide.runPython("import micropip")
  for (const pkg of packages) {
    self.postMessage({
      type: 'loading',
      msg: `Loading ${pkg}`,
      id: cell_id
    });
    await self.pyodide.runPythonAsync(`
      await micropip.install('${pkg}', keep_going=True);
    `);
  }
  console.log("Packages loaded!");
}

const autodetect_deps_code = `
import json
try:
    from panel.io.mime_render import find_requirements
except Exception:
    from panel.io.mime_render import find_imports as find_requirements
json.dumps(find_requirements(msg.to_py()['code']))`

const exec_code = `
from functools import partial
from panel.io.pyodide import pyrender
from js import console

msg = msg.to_py()
code = msg['code']
stdout_cb = partial(sendStdout, msg['id'])
stderr_cb = partial(sendStderr, msg['id'])
target = f"output-{msg['id']}"
pyrender(code, stdout_cb, stderr_cb, target)`

const onload_code = `
msg = msg.to_py()
if msg['mime'] == 'application/bokeh':
    from panel.io.pyodide import _link_docs_worker
    from panel.io.state import state
    doc = state.cache[f"output-{msg['id']}"]
    _link_docs_worker(doc, sendPatch, msg['id'], 'js')`

const patch_code = `
from panel import state

try:
    from pane.io.pyodide import _convert_json_patch
    patch = _convert_json_patch(msg.patch)
except:
    patch = msg.patch.to_py()
doc = state.cache[f"output-{msg.id}"]
doc.apply_json_patch(patch, setter='js')`

const MESSAGES = {
  patch: patch_code,
  execute: exec_code,
  rendered: onload_code
}

self.onmessage = async (event) => {
  let resolveExecution, rejectExecution;
   const executing = new Promise(function(resolve, reject){
    resolveExecution = resolve;
    rejectExecution = reject;
  });

  const prev_msg = QUEUE[0]
  const msg = {...event.data, executing}
  QUEUE.unshift(msg)

  if (prev_msg) {
    self.postMessage({
      type: 'loading',
      msg: 'Awaiting previous cells',
      id: msg.id
    });
    await prev_msg.executing
  }

  // Init pyodide
  if (self.pyodide == null) {
    self.postMessage({
      type: 'loading',
      msg: 'Loading pyodide',
      id: msg.id
    });
    await loadApplication(msg.id, msg.path)
    self.postMessage({
      type: 'loaded',
      id: msg.id
    });
  }

  // Handle message
  if (!MESSAGES.hasOwnProperty(msg.type)) {
    console.warn(`Service worker received unknown message type '${msg.type}'.`)
    resolveExecution()
    self.postMessage({
      type: 'idle',
      id: msg.id
    });
    return
  }


  try {
    self.pyodide.globals.set('msg', msg)
    let out = await self.pyodide.runPythonAsync(MESSAGES[msg.type])
    resolveExecution()
    if (out == null) {
      out = new Map()
    }
    if (out.has('content')) {
      self.postMessage({
        type: 'render',
        id: msg.id,
        content: out.get('content'),
        mime: out.get('mime_type')
      });
    }
    if (out.has('stdout') && out.get('stdout').length) {
      self.postMessage({
        type: 'stdout',
        content: out.get('stdout'),
        id: msg.id
      });
    }
    if (out.has('stderr') && out.get('stderr').length) {
      self.postMessage({
        type: 'stderr',
        content: out.get('stderr'),
        id: msg.id
      });
    }
    self.postMessage({
      type: 'idle',
      id: msg.id,
      uuid: msg.uuid
    });
  } catch (e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'error',
      traceback: traceback,
      msg: tblines[tblines.length-2],
      id: msg.id
    });
    resolveExecution()
    throw(e)
  }
}