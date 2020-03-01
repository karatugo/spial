import os
from flask import Flask, redirect, render_template, request, session, url_for
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = "./uploaded-fa-files"
app.config['PDB_UPLOAD_FOLDER'] = "./static/jsmol/data/"
head_title = "Spial"


@app.route('/', methods=["POST", "GET"])
def index():
    return render_template("index.html", head_title=head_title)

def is_fa_file(f):
    # TODO
    return True

def save_alignment_file_to_fs(f):
    filename = secure_filename(f.filename)
    f.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

@app.route('/spial', methods=["POST"])
def spial():
    alignment_a = request.files["alignment_a"]
    if alignment_a and is_fa_file(alignment_a.filename):
        save_alignment_file_to_fs(alignment_a)

    alignment_b = request.files["alignment_b"]
    if alignment_b and is_fa_file(alignment_b.filename):
        save_alignment_file_to_fs(alignment_b)
    
    consensus_threshold = request.form["consensus_threshold"]
    specificty_threshold = request.form["specificty_threshold"]

    return f"consensus: {consensus_threshold}, specificity: {specificty_threshold}"

if __name__ == '__main__':
    app.run(debug=True)

