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

    def sqlite_to_frame(self, sqlite_file):
        # where <path> is relative:
        engine = create_engine(f"sqlite:///{sqlite_file}")
        # todo do not leave this where filter in!!!
        # my_query = "SELECT * FROM parts_sequences where type = 'insertion'"
        my_query = "SELECT * FROM parts_sequences"
        results_frame = pd.read_sql_query(my_query, engine)
        engine.dispose()
        return results_frame

    def postgres_to_frame(self):
        # assumes that ssh tunnel is established
        # database_uri = f"postgresql+psycopg2://mam:{self.secrets_dict['db_pass']}@localhost:1111/felix"
        database_uri = f"postgresql+psycopg2://{self.secrets_dict['db_user']}@localhost:1111/felix"
        engine = create_engine(database_uri)
        my_query = "SELECT * FROM parts_sequences"
        results_frame = pd.read_sql_query(my_query, engine)
        engine.dispose()
        return results_frame
