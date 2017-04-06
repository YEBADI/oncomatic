from flask import Flask, render_template, url_for

# Define the WSGI application object
app = Flask(__name__)


# Configurations
app.config.from_pyfile('../config.py')

import views


