# not using ORM yet

# remove hard-coded secrets file path
# remove all hardcoded sql configuration


# from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Date
# from sqlalchemy.sql import text
# import sqlalchemy

# from sqlalchemy import create_engine
# from sqlalchemy import inspect
# from sqlalchemy.ext.declarative import declarative_base
# import pandas as pd
# import pprint

import pandas as pd
import yaml
from sqlalchemy import create_engine

secrets_file = '../local/secrets.yaml'

with open(secrets_file, 'r') as stream:
    try:
        secrets_dict = yaml.safe_load(stream)
        print(secrets_dict)
    except yaml.YAMLError as exc:
        print(exc)

# Base = declarative_base()

DATABASE_URI = f"postgresql+psycopg2://mam:{secrets_dict['dbpass']}@localhost:1111/felix"

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
