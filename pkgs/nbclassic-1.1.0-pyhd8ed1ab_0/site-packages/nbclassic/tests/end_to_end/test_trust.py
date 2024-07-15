from .utils import EDITOR_PAGE

trust_test_code = """
import IPython.display
IPython.display.Javascript("window.javascriptExecuted = true")
"""


def test_trust(notebook_frontend):
    def save_notebook():
        notebook_frontend.evaluate("() => Jupyter.notebook.save_notebook()", page=EDITOR_PAGE)

    # Add a cell that executes javascript
    notebook_frontend.add_cell(index=0, content=trust_test_code)
    notebook_frontend.execute_cell(0)
    notebook_frontend.wait_for_cell_output(0)
    # Make sure the JavaScript executed
    assert notebook_frontend.evaluate("() => window.javascriptExecuted", page=EDITOR_PAGE) == True
    # Check that we see the 'Trusted' text on the page
    trusted = notebook_frontend.locate("#notification_trusted", page=EDITOR_PAGE)
    assert trusted.get_inner_text() == "Trusted"
    save_notebook()

    # refresh the page, should be trusted
    notebook_frontend.reload(EDITOR_PAGE)
    notebook_frontend.wait_for_kernel_ready()
    trusted = notebook_frontend.locate("#notification_trusted", page=EDITOR_PAGE)
    assert trusted.get_inner_text() == "Trusted"
    notebook_frontend.wait_for_cell_output(0)
    assert notebook_frontend.evaluate("() => window.javascriptExecuted", page=EDITOR_PAGE) == True
