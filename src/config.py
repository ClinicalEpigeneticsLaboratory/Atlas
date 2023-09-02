import json
from os.path import exists


class ConfigFileNotFound(Exception):
    pass


def load_config(path: str) -> dict:
    if not exists(path):
        raise ConfigFileNotFound

    with open(path, "r") as file:
        config_file = json.load(file)
        return config_file


def export_config(config: dict, path: str) -> None:
    with open(path, "w") as file:
        file.write(json.dumps(config))


def add_filters(primary_site: str, diagnosis: str, config_dict: dict) -> dict:
    """

    :param primary_site:
    :param diagnosis:
    :param config_dict:
    :return: Dict
    """
    new_filters = [
        {
            "op": "=",
            "content": {
                "field": "cases.diagnoses.tissue_or_organ_of_origin",
                "value": primary_site,
            },
        },
        {
            "op": "=",
            "content": {
                "field": "cases.diagnoses.primary_diagnosis",
                "value": diagnosis,
            },
        },
    ]

    config_dict["FILTERS"]["content"].extend(new_filters)
    return config_dict
