#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import path as os_path
from flask import Flask, render_template, url_for, request, make_response
import argparse
import glob
import json
import io
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from fourgp_speclib import SpectrumLibrarySqlite

app = Flask(__name__)

# Set path to workspace where we expect to find libraries of spectra
our_path = os_path.split(os_path.abspath(__file__))[0]
workspace = os_path.join(our_path, "../../../4most-4gp-scripts/workspace")

# Read input parameters
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--library-path', required=False, dest='path', default=workspace,
                    help="Path to the collection of spectrum libraries we want to browse.")
parser.add_argument('--public',
                    required=False,
                    action='store_true',
                    dest="public",
                    help="Make this python/flask instance publicly visible on the network.")
parser.add_argument('--private',
                    required=False,
                    action='store_false',
                    dest="public",
                    help="Make this python/flask instance only visible on localhost (default).")
parser.set_defaults(public=False)
args = parser.parse_args()


# Index of all the libraries in this workspace
@app.route("/")
def library_index():
    # Fetch a list of all sub-directories inside this workspace -- each directory is a SpectrumLibrary
    libraries = glob.glob(os_path.join(args.path, "*"))
    libraries.sort()

    # For each library, look up how many spectra are inside it, and create a dictionary of properties
    library_info = []
    for item in libraries:
        if os_path.isdir(item):
            name = os_path.split(item)[1]
            x = SpectrumLibrarySqlite(path=item)
            library_info.append({
                'name': name,
                'url': url_for('library_search', library=name),
                'item_count': len(x)
            })
            x.close()
            del x

    # Render list of SpectrumLibraries into HTML
    return render_template('index.html', path=args.path, libraries=library_info)


# Search for spectra within a particular spectrum library
@app.route("/library/<library>", methods=("GET", "POST"))
def library_search(library):
    self_url = url_for("library_search", library=library)
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    metadata_keys = [str(i) for i in x._metadata_fields]
    metadata_keys.sort()

    # Fetch search constraints from POST data
    search = {"minima": {}, "maxima": {}}  # This structure contains strings which we put into the search form
    constraints = {}  # This structure contains the float / string constraints which we pass to the SpectrumLibrary
    # Loop over each of the metadata items we can be constrained on
    for item in metadata_keys:
        lower_limit = str(request.form.get("min_{}".format(item), ""))
        upper_limit = str(request.form.get("max_{}".format(item), ""))
        search['minima'][item] = lower_limit
        search['maxima'][item] = upper_limit
        # We only have a constraint on this metadata item if we have POST data which is not blank
        have_constraint = (lower_limit != "") or (upper_limit != "")

        # If we do have a constraint, test whether it's a string constraint, or a numeric constraint
        if have_constraint:
            lower_limit_float = 0
            upper_limit_float = 1e9
            # Test whether we get an exception when we try converting constraint to floats
            try:
                if lower_limit != "":
                    lower_limit_float = float(lower_limit)
                if upper_limit != "":
                    upper_limit_float = float(upper_limit)
                string_constraint = False
            except ValueError:
                string_constraint = True

            # Create a new numeric metadata constraint
            if not string_constraint:
                constraints[item] = ((lower_limit_float, upper_limit_float) if lower_limit_float != upper_limit_float
                                     else lower_limit_float)

            # Create a new string metadata constraint
            else:
                if lower_limit is None:
                    lower_limit = ""
                if upper_limit is None or upper_limit == "":
                    upper_limit = "zzzzzzzzz"
                constraints[item] = (lower_limit, upper_limit) if lower_limit != upper_limit else lower_limit

    # Search the SpectrumLibrary for matching spectra
    spectrum_ids = [i['specId'] for i in x.search(**constraints)]
    result_count = len(spectrum_ids)

    # Show a maximum of 100 results
    if len(spectrum_ids) > 100:
        spectrum_ids = spectrum_ids[:100]
    results = x.get_metadata(ids=spectrum_ids)

    # Add spectrum_id into each spectrum's metadata -- the HTML template needs this so we can link to spectrum viewer
    for i in range(len(spectrum_ids)):
        results[i]["spectrum_id"] = spectrum_ids[i]
    return render_template('library.html', path=args.path, library=library, metadata_keys=metadata_keys,
                           search=search, results=results, result_count=result_count, self_url=self_url)


# Display a particular spectrum
@app.route("/spectrum/<library>/<spec_id>", methods=("GET", "POST"))
def spectrum_view(library, spec_id):
    lambda_min = 3600
    lambda_max = 9600
    try:
        lambda_min = float(request.form.get("lambda_min"))
    except (TypeError, ValueError):
        pass
    try:
        lambda_max = float(request.form.get("lambda_max"))
    except (TypeError, ValueError):
        pass

    parent_url = url_for("library_search", library=library)
    self_url = url_for("spectrum_view", library=library, spec_id=spec_id)
    txt_url = url_for("spectrum_txt", library=library, spec_id=spec_id)
    data_url = url_for("spectrum_json", library=library, spec_id=spec_id)
    png_url = url_for("spectrum_png", library=library, spec_id=spec_id, lambda_min=lambda_min, lambda_max=lambda_max)

    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)

    metadata_keys = x.list_metadata_fields()
    metadata_keys.sort()
    metadata = x.get_metadata(ids=int(spec_id))[0]
    metadata["spectrum_id"] = spec_id

    return render_template('spectrum.html', path=args.path, library=library, metadata_keys=metadata_keys,
                           parent_url=parent_url, metadata=metadata,
                           txt_url=txt_url, data_url=data_url, png_url=png_url,
                           self_url=self_url, lambda_min=lambda_min, lambda_max=lambda_max)


# Output a particular spectrum as a JSON file
@app.route("/spectrum_json/<library>/<spec_id>")
def spectrum_json(library, spec_id):
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    spectrum = x.open(ids=int(spec_id)).extract_item(0)

    data = list(zip(spectrum.wavelengths, spectrum.values))
    return json.dumps(data)


# Output a particular spectrum as a JSON file
@app.route("/spectrum_txt/<library>/<spec_id>")
def spectrum_txt(library, spec_id):
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    spectrum = x.open(ids=int(spec_id)).extract_item(0)
    data = np.asarray(list(zip(spectrum.wavelengths, spectrum.values, spectrum.value_errors)))

    txt_output = io.StringIO()
    np.savetxt(txt_output, data)
    response = make_response(txt_output.getvalue())
    response.headers['Content-Type'] = 'text/plain'
    return response


# Output a particular spectrum as a png file
@app.route("/spectrum_png/<library>/<spec_id>/<lambda_min>/<lambda_max>")
def spectrum_png(library, spec_id, lambda_min, lambda_max):
    path = os_path.join(args.path, library)
    x = SpectrumLibrarySqlite(path=path)
    spectrum = x.open(ids=int(spec_id)).extract_item(0)

    fig = Figure(figsize=(16, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Wavelength / A')
    ax.set_ylabel('Value')
    ax.set_xlim([float(lambda_min), float(lambda_max)])
    ax.grid(True)
    ax.plot(spectrum.wavelengths, spectrum.values)
    canvas = FigureCanvas(fig)
    png_output = io.StringIO()
    canvas.print_png(png_output)
    response = make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


if __name__ == "__main__":
    app.run(host="0.0.0.0" if args.public else "127.0.0.1")
