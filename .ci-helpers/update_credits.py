import pathlib
import re
import textwrap
import warnings

import requests
from rst_include import rst_inc


def generate_zenodo():
    """Generates `zenodo.rst` file with BibTeX citation
    Adapted from: https://astrodata.nyc/posts/2021-04-23-zenodo-sphinx/"""

    CONCEPT_DOI = "592480"  # See: https://help.zenodo.org/#versioning
    zenodo_path = pathlib.Path("docs/resources/zenodo.rst")

    try:
        headers = {"accept": "application/x-bibtex"}
        response = requests.get(
            f"https://zenodo.org/api/records/{CONCEPT_DOI}", headers=headers
        )
        response.encoding = "utf-8"
        citation = re.findall("@software{(.*)\,", response.text)
        zenodo_record = (
            f".. |ZENODO| replace:: {citation[0]}\n\n"
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


if __name__ == "__main__":
    generate_zenodo()
    rst_inc(source='docs/resources/credits_template.rst', target='docs/resources/credits.rst')
    rst_inc(source='README_TEMPLATE.rst', target='README.rst')
