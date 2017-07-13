#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path as os_path
from flask import Flask, render_template, url_for, request
import argparse
import glob
import json

from fourgp_speclib import SpectrumLibrarySqlite

app = Flask(__name__)

# Set path to workspace where we expect to find libraries of spectra
our_path = os_path.split(os_path.abspath(__file__))[0]
workspace = os_path.join(our_path, "../../../4most-4gp-scripts/workspace")

# Read input parameters
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--library-path', required=False, dest='path', default=workspace,
                    help="Path to the collection of spectrum libraries we want to browse.")
args = parser.parse_args()


# Index of all the libraries in this workspace
@app.route("/")
def library_index():
    libraries = glob.glob(os_path.join(args.path, "*"))
    libraries.sort()

    library_info = []
    for item in libraries:
        name = os_path.split(item)[1]
        x = SpectrumLibrarySqlite(path=item)
        library_info.append({
            'name': name,
            'url': url_for('library_search', library=name),
            'item_count': len(x.search())
        })
        x.close()
        del x
    return render_template('index.html', path=args.path, libraries=library_info)


# Search a particular spectrum library
@app.route("/library/<library>")
def library_search(library):
    self_url = url_for("library_search", library=library)
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    metadata_keys = x._metadata_fields
    metadata_keys.sort()
    search = {"minima": {}, "maxima": {}}
    constraints = {}
    for item in metadata_keys:
        search['minima'][item] = ""
        search['maxima'][item] = ""
    spectrum_ids = [i['specId'] for i in x.search(**constraints)]
    result_count = len(spectrum_ids)
    if len(spectrum_ids) > 100:
        spectrum_ids = spectrum_ids[:100]
    results = x.get_metadata(ids=spectrum_ids)
    for i in range(len(spectrum_ids)):
        results[i]["spectrum_id"] = spectrum_ids[i]
    return render_template('library.html', path=args.path, library=library, metadata_keys=metadata_keys,
                           search=search, results=results, result_count=result_count, self_url=self_url)


# Search a particular spectrum
@app.route("/spectrum/<library>/<spec_id>")
def spectrum_view(library, spec_id):
    parent_url = url_for("library_search", library=library)
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    metadata_keys = x._metadata_fields
    metadata_keys.sort()
    metadata = x.get_metadata(ids=int(spec_id))[0]
    metadata["spectrum_id"] = spec_id
    return render_template('spectrum.html', path=args.path, library=library, metadata_keys=metadata_keys,
                           parent_url=parent_url, metadata=metadata )


# Output a particular spectrum as a JSON file
@app.route("/spectrum_json/<library>/<spec_id>")
def spectrum_json(library, spec_id):
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    spectrum = x.open(ids=int(spec_id)).extract_item(0)
    data = zip(spectrum.wavelengths, spectrum.values)
    return json.dumps(data)


if __name__ == "__main__":
    app.run()
