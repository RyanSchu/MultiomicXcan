import argparse as ap
import sqlite3

parser = ap.ArgumentParser()
parser.add_argument('--db', help = 'file path of the target db file to be parsed')
parser.add_argument('--out', help='File out with path')
args = parser.parse_args()


conn = sqlite3.connect(args.db)
cursor = conn.execute("SELECT * FROM weights ORDER BY gene;")
with open(args.out, "a") as file:
        for row in cursor:
                file.write('\t'.join(str(e) for e in row) + '\n')
conn.close()
