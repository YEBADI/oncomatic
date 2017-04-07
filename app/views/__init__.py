from flask import render_template, request
from app import app 
import subprocess

@app.route('/', methods=['GET'])
def home_page():
    return render_template('home.html')

@app.route('/about', methods=['GET'])
def about_page():
    return render_template('about.html')




@app.route('/result', methods=['POST'])
def results_page():
    print(request.form)
    tumortype = request.form['tumor_type']
    pipeline = request.form['pipeline']
    genenumber = request.form['genenumber']
    pickgenes = request.form['pickgenes']
    genelist = request.form['genelist']

    if tumortype == 'BRCA,subtyping':
      subprocess.call(['Rscript', 'app/scripts/BRCA_oncoprint', 
                   'yes', pipeline, genenumber, pickgenes, genelist])
      
      return render_template('results_subtyped.html', 
                         pickgenes=pickgenes, genelist=genelist, 
                         genenumber = genenumber, tumortype = tumortype, 
                         pipeline = pipeline)

    else:
      subprocess.call(['Rscript', 'app/scripts/all_tumors_oncoprint', 
                   tumortype, pipeline, genenumber, pickgenes, genelist])

      return render_template('results.html', 
                           pickgenes=pickgenes, genelist=genelist, 
                           genenumber = genenumber, tumortype = tumortype, 
                           pipeline = pipeline)

# Sample HTTP error handling
@app.errorhandler(404)
def not_found(error):
    return render_template('404.html'), 404
