# -*- coding: utf-8 -*-

# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import os
from contextlib import ContextDecorator

import pymysql



class ProteinDatabaseCursor(ContextDecorator):
    def __init__(self):
        self.connection = None

    def __enter__(self):
        config_file = os.getenv("METAGENOMICS_DB_CONFIG")
        if config_file is None:
            raise ValueError(
                "Please set the value of the env variable METAGENOMICS_DB_CONFIG"
            )

        with open(config_file) as creds_file:
            configuration = json.load(creds_file)

        self.connection = pymysql.connect(
            host=configuration.get("host"),
            user=configuration.get("user"),
            password=configuration.get("password"),
            port=configuration.get("port"),
            cursorclass=pymysql.cursors.DictCursor,
        )
        cursor = self.connection.cursor()
        return cursor

    def __exit__(self, *exc):
        self.connection.close()
        return False
