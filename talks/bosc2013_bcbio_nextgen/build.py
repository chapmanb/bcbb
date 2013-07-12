"""Convert IPython notebook into HTML reveal.js slides.

Overrides some standard functionality:
  - Uses a local reveal.js with customized CSS.
  - Avoids writing out the standard IPython reveal CSS overrides in favor of the default style.
"""
from IPython.nbconvert.exporters import RevealExporter
from IPython.config import Config

from IPython.nbformat import current as nbformat

infile = "chapmanb_bosc2013_bcbio.ipynb"
outfile = "chapmanb_bosc2013_bcbio.html"

notebook = open(infile).read()
notebook_json = nbformat.reads_json(notebook)

c = Config({'RevealHelpTransformer': {'enabled': True,
                                      'url_prefix':'../reveal.js',},
            "CSSHTMLHeaderTransformer": {'enabled': False}
            })

exportHtml = RevealExporter(config=c)
(body,resources) = exportHtml.from_notebook_node(notebook_json)

with open(outfile, "w") as out_handle:
    in_css_override = False
    for line in body.encode('utf-8').split("\n"):
        if line.startswith("/* Overrides of notebook CSS"):
            in_css_override = True
        if in_css_override:
            if line.startswith("</style>"):
                in_css_override = False
        if not in_css_override:
            out_handle.write(line + "\n")
