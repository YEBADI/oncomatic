from flask import render_template, request
from app import app 

@app.route('/', methods=['GET'])
def home_page():
    return render_template('home.html')

@app.route('/about', methods=['GET'])
def about_page():
    return render_template('about.html')

@app.route('/result', methods=['POST'])
def results_page():
    print(request.form)
    name = request.form['first_name']
    telephone = request.form['telephone']
    return render_template('results.html', 
                           name = name, tel = telephone)


# Sample HTTP error handling
@app.errorhandler(404)
def not_found(error):
    return render_template('404.html'), 404
