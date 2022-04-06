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
