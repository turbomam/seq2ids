# not using ORM yet

# get that pw out of there


# from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Date
# from sqlalchemy.sql import text
# import sqlalchemy

# from sqlalchemy import create_engine
# from sqlalchemy import inspect
# from sqlalchemy.ext.declarative import declarative_base
# import pandas as pd
# import pprint

import pandas as pd
from sqlalchemy import create_engine

# Base = declarative_base()

DATABASE_URI = 'postgresql+psycopg2://mam:PASSWORD@localhost:1111/felix'

engine = create_engine(DATABASE_URI)

# inspector = inspect(engine)
# temp = inspector.get_columns('parts_sequences')
# pprint.pprint(temp)

myQuery = "SELECT id, file FROM parts_sequences limit 9"

# with engine.connect() as con:
#     rs = con.execute(myQuery)
#
#     for row in rs:
#         print(row)

df = pd.read_sql_query(myQuery, engine)

print(df)

# Base.metadata.drop_all(engine)

engine.dispose()

# id
# part_id
# file
# seq_name
# type
# date_added
# sequence
