from flask import Flask, render_template
from pdb import get_first_pdb_id

app = Flask(__name__)

@app.route('/')
def index():
    first_pdb_id = get_first_pdb_id()
    return render_template('index.html', first_pdb_id=first_pdb_id)

if __name__ == '__main__':
    app.run(debug=True)







