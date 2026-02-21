"""Script for updating `credits.rst` and `README.rst` between releases."""

import pathlib
import re
import textwrap
import warnings
from datetime import date

import requests

INCLUDE_PATTERN = re.compile(r"^\.\. include:: (.+)$", re.MULTILINE)


def resolve_rst_includes(source, target):
    """Resolve ``.. include::`` directives by inlining referenced files.

    Parameters
    ----------
    source : str
        Path to the RST template containing include directives.
    target : str
        Path to write the resolved output.
    """
    source_path = pathlib.Path(source)
    text = source_path.read_text(encoding="utf-8")

    def _replace(match):
        include_path = source_path.parent / match.group(1)
        return include_path.read_text(encoding="utf-8")

    resolved = INCLUDE_PATTERN.sub(_replace, text)
    pathlib.Path(target).write_text(resolved, encoding="utf-8")


def generate_zenodo():
    """Generates `zenodo.rst` file with BibTeX citation
    Adapted from: https://astrodata.nyc/posts/2021-04-23-zenodo-sphinx/"""

    CONCEPT_DOI = "592480"  # See: https://help.zenodo.org/#versioning
    zenodo_path = pathlib.Path("docs/resources/zenodo.rst")
    year = date.today().year

    try:
        headers = {"accept": "application/x-bibtex"}
        response = requests.get(
            f"https://zenodo.org/api/records/{CONCEPT_DOI}", headers=headers
        )
        response.encoding = "utf-8"
        citation = re.findall(r"@software{(.*)\,", response.text)
        specific_doi = citation[0].lstrip(f"kerzendorf_wolfgang_{year}_")
        doi_org_url = f"https://doi.org/10.5281/zenodo.{specific_doi}"
        doi_badge = f"https://img.shields.io/badge/DOI-10.5281/zenodo.{specific_doi}-blue"
        zenodo_record = (
            f".. |CITATION| replace:: {citation[0]}\n\n"
            f".. |DOI_BADGE| image:: {doi_badge}\n"
            f"                 :target: {doi_org_url}\n\n"
            ".. code-block:: bibtex\n\n"
            + textwrap.indent(response.text, " " * 4)
        )

    except Exception as e:
        warnings.warn(
            "Failed to retrieve Zenodo record for TARDIS: " f"{str(e)}"
        )

        not_found_msg = """
                        Couldn"t retrieve the TARDIS software citation from Zenodo. Get it 
                        directly from `this link <https://zenodo.org/record/{CONCEPT_DOI}>`_    .
                        """

        zenodo_record = (
            ".. |ZENODO| replace:: <TARDIS SOFTWARE CITATION HERE> \n\n"
            ".. warning:: \n\n" + textwrap.indent(not_found_msg, " " * 4)
        )

    with open(zenodo_path, "w") as f:
        f.write(zenodo_record)

    print(zenodo_record)


def main():
    generate_zenodo()

    resolve_rst_includes(
        source="docs/resources/credits_template.rst",
        target="docs/resources/credits.rst",
    )

    resolve_rst_includes(
        source="README_TEMPLATE.rst",
        target="README.rst",
    )


if __name__ == "__main__":
    main()
