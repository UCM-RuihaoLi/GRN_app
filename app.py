from flask import Flask, render_template, Response, request
import subprocess
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import random
import math
import copy
from scipy.integrate import solve_ivp
from distance_functions import *
from mutation_functions import *
from dynamics import *
from utility_functions import *
from population_functions import *
from GRN_Expanded_Combinatorial import GRN
from Combine_Redundant_Attractors import *
from flask_sslify import SSLify

app = Flask(__name__)
sslify = SSLify(app)
input_value = None  # Store the input value globally

@app.route('/')
def index():
    return render_template('GRN_Simulator_test.html')

@app.route('/set-input', methods=['POST'])
def set_input():
    global input_value
    input_value = request.form['input_value']
    return '', 204  # No Content

@app.route('/stream-data')
def stream_data():
    global input_value

    def generate():
        with subprocess.Popen(['python', 'GRN_simulator_test.py', input_value], stdout=subprocess.PIPE, text=True) as process:
            for line in process.stdout:
                if line.strip() == 'Process completed.':
                    yield f"data: {line}\n\n"  # Signal completion
                    break
                yield f"data: {line}\n\n"
            process.stdout.close()
            process.wait()
        yield "data: Attractor calculation completed.\n\n"  # Signal completion

    return Response(generate(), mimetype='text/event-stream')

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
