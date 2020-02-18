from flask import Flask, redirect, render_template, request, session, url_for

app = Flask(__name__)
head_title = "Spial"


@app.route('/', methods=["POST", "GET"])
def index():
    return render_template("index.html", head_title=head_title)


if __name__ == '__main__':
    app.run(debug=True)

