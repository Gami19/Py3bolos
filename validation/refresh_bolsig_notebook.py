"""Execute the BOLSIG+ notebook, optionally exporting its numerical tables.

Run ``python validation/refresh_bolsig_notebook.py --export`` from the source
checkout. The notebook's commented-out export cell is preserved. With --export,
that cell is executed temporarily after the calculation, then removed before
saving the notebook. This provides explicit export control without changing
the user's notebook setting.
"""

import argparse
from pathlib import Path
import sys

from jupyter_client import KernelManager
import nbformat
from nbclient import NotebookClient


def main():
    """Run the notebook with this Python environment and save its outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = root / "validate_bolsig.ipynb"
    notebook = nbformat.read(path, as_version=4)
    temporary_export = False
    if args.export:
        cell = notebook.cells[-1]
        if cell.cell_type != "code" or not cell.source.startswith("# tables ="):
            raise ValueError("Expected the existing commented export cell.")
        source = "\n".join(
            line[2:] if line.startswith("# ") else line
            for line in cell.source.splitlines()
        )
        notebook.cells.append(nbformat.v4.new_code_cell(source))
        temporary_export = True
    manager = KernelManager(kernel_name="python3")
    manager.kernel_spec.argv = [
        sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"
    ]
    client = NotebookClient(
        notebook, km=manager, timeout=300,
        resources={"metadata": {"path": str(root)}},
        on_cell_executed=lambda cell_index, **kwargs: print(
            f"Executed cell {cell_index + 1}", flush=True
        ),
    )
    completed = False
    try:
        client.execute()
        completed = True
    finally:
        if temporary_export:
            notebook.cells.pop()
        if manager.has_kernel:
            manager.shutdown_kernel(now=True)
        notebook.metadata["bolos_validation"] = {
            "local_exponential_sign": "corrected",
            "execution_completed": completed,
            "tables_exported": completed and args.export,
        }
        nbformat.validate(notebook)
        nbformat.write(notebook, path)
    print("Notebook executed; CSV/JSON export:", args.export)


if __name__ == "__main__":
    main()
