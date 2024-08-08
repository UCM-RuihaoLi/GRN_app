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
import paramiko
import getpass
from dotenv import load_dotenv

app = Flask(__name__)
sslify = SSLify(app)
input_value = None  # Store the input value globally

# Load environment variables from .env file
load_dotenv()

def run_computation_on_vm(data):
    host = "128.84.29.86"
    port = 22
    username = "ubuntu"
    private_key_path = os.getenv("PRIVATE_KEY_PATH")
    passphrase = os.getenv("PRIVATE_KEY_PASSPHRASE")

    # Prepare data to send
    data_str = ' '.join(map(str, data))

    # Set up the SSH client
    ssh_client = paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # Load private key with passphrase
    private_key = paramiko.RSAKey(filename=private_key_path, password=passphrase)

    ssh_client.connect(hostname=host, port=port, username=username, pkey=private_key)

    # Activate conda environment and run the computation script
    command = f"source ~/miniconda3/bin/activate && conda activate base && echo '{data_str}' | python ~/test.py"
    stdin, stdout, stderr = ssh_client.exec_command(command)
    result = stdout.read().decode().strip()
    error = stderr.read().decode().strip()
    ssh_client.close()

    if error:
        print(f"Error: {error}")
    return result

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
