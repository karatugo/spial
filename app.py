import flask
import os
from dotenv import load_dotenv
from flask import Flask, redirect, render_template, request, session, url_for
from flask_sqlalchemy import SQLAlchemy
from werkzeug.utils import secure_filename
from src.spial import get_aa_conservation, get_sdp

app = Flask(__name__)
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///test.db"
db = SQLAlchemy(app)

load_dotenv()
app.secret_key = os.getenv("SECRET_KEY")
app.config['UPLOAD_FOLDER'] = "./uploaded-fa-files"
app.config['PDB_UPLOAD_FOLDER'] = "./static/jsmol/data/"

FA_EXTENSIONS = {"fa"}
PDB_EXTENSIONS = {"pdb"}
head_title = "Spial"

class Calculation(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    conservation_dict1 = db.Column(db.String(200), nullable=False)
    conservation_dict2 = db.Column(db.String(200), nullable=False)
    consensus_threshold = db.Column(db.DECIMAL, nullable=False)
    specificty_threshold = db.Column(db.DECIMAL, nullable=False)
    sdp_dict = db.Column(db.Text, nullable=False)
    pdb_filename = db.Column(db.String(200), nullable=False)

    def __repr__(self):
        return "<Calculation {}>".format(self.id)


@app.route('/', methods=["POST", "GET"])
def index():
    pdb_filename = None
    if 'pdb_filename' in flask.session:
        pdb_filename = flask.session['pdb_filename']

    return render_template("index.html", 
                           head_title=head_title,
                           pdb_filename=pdb_filename,
                           sdp_dict=None)

def is_fa_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in FA_EXTENSIONS

def is_pdb_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in PDB_EXTENSIONS

def save_file_to_fs(upload_folder, f):
    filename = secure_filename(f.filename)
    f.save(os.path.join(upload_folder, filename))

@app.route('/spial', methods=["POST"])
def spial():
    session['key'] = 'value'
    alignment_a = request.files["alignment_a"]
    if alignment_a and is_fa_file(alignment_a.filename):
        save_file_to_fs(app.config['UPLOAD_FOLDER'], alignment_a)

    alignment_b = request.files["alignment_b"]
    if alignment_b and is_fa_file(alignment_b.filename):
        save_file_to_fs(app.config['UPLOAD_FOLDER'], alignment_b)
    
    pdb = request.files["pdb"]
    if pdb and is_pdb_file(pdb.filename):
        save_file_to_fs(app.config['PDB_UPLOAD_FOLDER'], pdb)
    
    consensus_threshold = request.form["consensus_threshold"]
    specificty_threshold = request.form["specificty_threshold"]

    conservation_dict1 = get_aa_conservation(alignment_a.filename)
    conservation_dict2 = get_aa_conservation(alignment_b.filename)
    sdp_dict = get_sdp(conservation_dict1,
                       conservation_dict2,
                       consensus_threshold,
                       specificty_threshold)

    new_calc = Calculation(conservation_dict1=str(conservation_dict1),
                           conservation_dict2=str(conservation_dict2),
                           consensus_threshold=consensus_threshold,
                           specificty_threshold=specificty_threshold,
                           sdp_dict=str(sdp_dict),
                           pdb_filename=pdb.filename)
    try:
        db.session.add(new_calc)
        db.session.commit()
        return redirect("/")
    except Exception as e:
        return f"There was an issue adding your task: {e}"


    flask.session["pdb_filename"] = pdb.filename
    return flask.redirect("/")

if __name__ == '__main__':
    app.run(debug=True)
