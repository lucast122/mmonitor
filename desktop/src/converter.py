import sqlite3

import pandas as pd
from mysql.connector import connect

with sqlite3.connect('mmonitor.sqlite3') as sqlite_con:
    tax_df = pd.read_sql_query("SELECT * FROM mmonitor", sqlite_con)
    meta_df = pd.read_sql_query("SELECT * FROM metadata", sqlite_con)

with connect(user='admin', password='admin', host='localhost', database='mmonitor') as mysql_con:
    cur = mysql_con.cursor()
    for values in tax_df.itertuples(index=False, name=None):
        cur.execute(f"INSERT INTO mmonitor VALUES {values}")
    for values in meta_df.itertuples(index=False, name=None):
        cur.execute(f"INSERT INTO metadata VALUES {values}")
    mysql_con.commit()
