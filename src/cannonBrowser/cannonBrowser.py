#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from os import path as os_path
from flask import Flask, render_template, url_for, request, make_response
import argparse
import datetime
import re
import glob
import json
import pickle
import StringIO
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

app = Flask(__name__)

# Set path to workspace where we expect to find Cannon output files
our_path = os_path.split(os_path.abspath(__file__))[0]
workspace = os_path.join(our_path, "../../../4most-4gp-scripts/output_data/cannon")

# Read input parameters
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-path', required=False, dest='path', default=workspace,
                    help="Path to the collection of Cannon output files we want to browse.")
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


# Index of all the Cannon runs in this workspace
@app.route("/")
def cannon_index():
    # Fetch a list of all the Cannon runs inside this workspace
    cannon_runs = glob.glob(os_path.join(args.path, "*.json"))
    cannon_runs.sort()

    # For each library, look up how many spectra are inside it, and create a dictionary of properties
    cannon_info = []
    for item in cannon_runs:
            name = os_path.split(item)[1]
            file_time = os_path.getmtime(item)
            cannon_info.append({
                'name': name,
                'url': url_for('cannon_run', cannon=name),
                'item_time': datetime.datetime.fromtimestamp(file_time).strftime('%Y-%m-%d %H:%M:%S')
            })

    # Render list of SpectrumLibraries into HTML
    return render_template('index.html', path=args.path, cannon_runs=cannon_info)


# Display summary of a particular Cannon run
@app.route("/cannon/<cannon>", methods=("GET", "POST"))
def cannon_run(cannon):
    self_url = url_for("cannon_run", cannon=cannon)
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)
    x = json.load(fp=open(path_json))

    y = pickle.load(file=open(path_cannon))
    vectorizer_terms = enumerate(y['vectorizer'].get_human_readable_label_vector().split(' + '))

    metadata = [
        {
            "key": "Training library",
            "value": x["train_library"]
         }, {
            "key": "Test library",
            "value": x["test_library"]
        }, {
            "key": "Description",
            "value": x["description"]
        }, {
            "key": "Labels",
            "value": str([str(i) for i in x["labels"]])
        }, {
            "key": "Start time",
            "value": datetime.datetime.fromtimestamp(x["start_time"]).strftime('%Y-%m-%d %H:%M:%S')
        }, {
            "key": "End time",
            "value": datetime.datetime.fromtimestamp(x["end_time"]).strftime('%Y-%m-%d %H:%M:%S')
        }, {
            "key": "Training time",
            "value": "{:.1f} sec".format(x["training_time"])
        }, {
            "key": "Test time (per spectrum)",
            "value": "{:.1f} sec".format((x["end_time"] - x["start_time"] - x["training_time"]) / len(x["stars"]))
        }, {
            "key": "Test time (total)",
            "value": "{:.1f} sec".format(x["end_time"] - x["start_time"] - x["training_time"])
        }, {
            "key": "Spectrum count",
            "value": len(x["stars"])
        }, {
            "key": "Wavelength range",
            "value": "{:d}A to {:d}A".format(int(x["wavelength_raster"][0]), int(x["wavelength_raster"][-1]))
        }, {
            "key": "Pixel count",
            "value": "{:d}".format(len(x["wavelength_raster"]))
        }, {
            "key": "Line list",
            "value": x["line_list"]
        }, {
            "key": "Assume scaled solar",
            "value": x["assume_scaled_solar"]
        }, {
            "key": "Tolerance",
            "value": x["tolerance"]
        }
    ]

    return render_template('cannon_run.html', path=args.path, cannon=cannon, metadata=metadata, self_url=self_url,
                           vectorizer_terms=vectorizer_terms)


# Display a particular coefficient spectrum
@app.route("/coefficient_spectrum/<cannon>/<term>", methods=("GET", "POST"))
def coefficient_spectrum(cannon, term):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))
    vectorizer_terms = y['vectorizer'].get_human_readable_label_vector().split(' + ')
    term_name = vectorizer_terms[int(term)]

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

    parent_url = url_for("cannon_run", cannon=cannon)
    self_url = url_for("coefficient_spectrum", cannon=cannon, term=term)
    txt_url = url_for("coefficient_spectrum_txt", cannon=cannon, term=term)
    data_url = url_for("coefficient_spectrum_json", cannon=cannon, term=term)
    png_url = url_for("coefficient_spectrum_png", cannon=cannon, term=term, lambda_min=lambda_min, lambda_max=lambda_max)

    return render_template('coefficient_spectrum.html', path=args.path, cannon=cannon, term_name=term_name,
                           parent_url=parent_url, txt_url=txt_url, data_url=data_url, png_url=png_url,
                           self_url=self_url, lambda_min=lambda_min, lambda_max=lambda_max)


# Output a particular coefficient spectrum as a JSON file
@app.route("/coefficient_spectrum_json/<cannon>/<term>")
def coefficient_spectrum_json(cannon, term):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))
    data = zip(y['dispersion'], y['theta'][:,int(term)])
    return json.dumps(data)


# Output a particular coefficient spectrum as a JSON file
@app.route("/coefficient_spectrum_txt/<cannon>/<term>")
def coefficient_spectrum_txt(cannon, term):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))
    data = zip(y['dispersion'], y['theta'][:,int(term)])

    txt_output = StringIO.StringIO()
    np.savetxt(txt_output, data)
    response = make_response(txt_output.getvalue())
    response.headers['Content-Type'] = 'text/plain'
    return response


# Output a particular coefficient spectrum as a png file
@app.route("/coefficient_spectrum_png/<cannon>/<term>/<lambda_min>/<lambda_max>")
def coefficient_spectrum_png(cannon, term, lambda_min, lambda_max):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))

    fig = Figure(figsize=(16, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Wavelength / A')
    ax.set_ylabel('Value')
    ax.set_xlim([float(lambda_min), float(lambda_max)])
    ax.grid(True)
    ax.plot(y['dispersion'], y['theta'][:,int(term)])
    canvas = FigureCanvas(fig)
    png_output = StringIO.StringIO()
    canvas.print_png(png_output)
    response = make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


# Display a particular scatter spectrum
@app.route("/scatter_spectrum/<cannon>", methods=("GET", "POST"))
def scatter_spectrum(cannon):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

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

    parent_url = url_for("cannon_run", cannon=cannon)
    self_url = url_for("scatter_spectrum", cannon=cannon)
    txt_url = url_for("scatter_spectrum_txt", cannon=cannon)
    data_url = url_for("scatter_spectrum_json", cannon=cannon)
    png_url = url_for("scatter_spectrum_png", cannon=cannon, lambda_min=lambda_min, lambda_max=lambda_max)

    return render_template('scatter_spectrum.html', path=args.path, cannon=cannon,
                           parent_url=parent_url, txt_url=txt_url, data_url=data_url, png_url=png_url,
                           self_url=self_url, lambda_min=lambda_min, lambda_max=lambda_max)


# Output a particular scatter spectrum as a JSON file
@app.route("/scatter_spectrum_json/<cannon>")
def scatter_spectrum_json(cannon):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))
    data = zip(y['dispersion'], y['s2'])
    return json.dumps(data)


# Output a particular scatter spectrum as a JSON file
@app.route("/scatter_spectrum_txt/<cannon>")
def scatter_spectrum_txt(cannon):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))
    data = zip(y['dispersion'], y['s2'])

    txt_output = StringIO.StringIO()
    np.savetxt(txt_output, data)
    response = make_response(txt_output.getvalue())
    response.headers['Content-Type'] = 'text/plain'
    return response


# Output a particular scatter spectrum as a png file
@app.route("/scatter_spectrum_png/<cannon>/<lambda_min>/<lambda_max>")
def scatter_spectrum_png(cannon, lambda_min, lambda_max):
    path_json = os_path.join(args.path, cannon)
    path_cannon = re.sub(".json", ".cannon", path_json)

    y = pickle.load(file=open(path_cannon))

    fig = Figure(figsize=(16, 6))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Wavelength / A')
    ax.set_ylabel('Value')
    ax.set_xlim([float(lambda_min), float(lambda_max)])
    ax.grid(True)
    ax.plot(y['dispersion'], y['s2'])
    canvas = FigureCanvas(fig)
    png_output = StringIO.StringIO()
    canvas.print_png(png_output)
    response = make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


if __name__ == "__main__":
    app.run(host="0.0.0.0" if args.public else "127.0.0.1")
