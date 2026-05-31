"""
Page generators for TARDIS documentation
Automatically generates index pages for tutorials and workflows.
"""
from pathlib import Path


def generate_tutorials_page(app):
    """Generate tutorials page with all tutorial notebooks."""
    notebooks = ""
    tutorials_path = Path("tutorials")

    for notebook in tutorials_path.rglob("*.ipynb"):
        if "tutorial" in notebook.name and "checkpoint" not in notebook.name:
            notebooks += f"\n* :doc:`{notebook.parent}/{notebook.stem}`"

    title = "Tutorials\n*********\n"
    description = "The following pages contain the TARDIS tutorials:\n\n"

    with open("tutorials.rst", mode="wt", encoding="utf-8") as f:
        f.write(f"{title}\n{description}\n{notebooks}")


def generate_workflows_page(app):
    """Generate workflows page with all workflow notebooks."""
    notebooks = ""
    workflows_path = Path("workflows")

    for notebook in workflows_path.rglob("*.ipynb"):
        if "workflow" in notebook.name and "checkpoint" not in notebook.name:
            notebooks += f"\n* :doc:`{notebook.parent}/{notebook.stem}`"

    title = "Workflows\n*********\n"
    description = "The following pages contain the TARDIS workflows:\n\n These examples are intended to help users explore specific modules within TARDIS, with the goal of supporting their individual scientific objectives."

    with open("workflows.rst", mode="wt", encoding="utf-8") as f:
        f.write(f"{title}\n{description}\n{notebooks}")


def setup(app):
    """Setup the page generators extension."""
    app.connect("builder-inited", generate_tutorials_page)
    app.connect("builder-inited", generate_workflows_page)
    
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
