# todo not using ORM yet
# todo remove all hardcoded sql configuration

from typing import Dict, Any

import pandas as pd
import yaml
from sqlalchemy import create_engine


class SeqsFromDb:
    def __init__(self):
        self.secrets_dict = Dict[str, Any]

    def get_secrets_dict(self, secrets_file):
        with open(secrets_file, 'r') as stream:
            try:
                self.secrets_dict = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

    def query_to_frame(self):
        database_uri = f"postgresql+psycopg2://mam:{self.secrets_dict['db_pass']}@localhost:1111/felix"
        engine = create_engine(database_uri)
        my_query = "SELECT * FROM parts_sequences"
        results_frame = pd.read_sql_query(my_query, engine)
        engine.dispose()
        return results_frame

# ---

# from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Date
# from sqlalchemy.sql import text
# import sqlalchemy

# from sqlalchemy import create_engine
# from sqlalchemy import inspect
# from sqlalchemy.ext.declarative import declarative_base
# import pandas as pd
# import pprint

# Base = declarative_base()

# inspector = inspect(engine)
# temp = inspector.get_columns('parts_sequences')
# pprint.pprint(temp)

# with engine.connect() as con:
#     rs = con.execute(myQuery)
#
#     for row in rs:
#         print(row)

# Base.metadata.drop_all(engine)

# id
# part_id
# file
# seq_name
# type
# date_added
# sequence
