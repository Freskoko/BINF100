import json


def load_config():
    with open("/home/henrik/school/BINF100/excerices/dna_programming/config.json") as file:
        return json.load(file)


CONFIG = load_config()