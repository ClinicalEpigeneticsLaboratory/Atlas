import json
import requests
from io import StringIO

import pandas as pd


def request_gdc_service(
    fields: list,
    filters: dict,
    n_records: int = 10000,
) -> pd.DataFrame:
    """
    Function requests GDC service to get files described in config file.

    :param fields:
    :param filters:
    :param n_records:
    :return: pd.DataFrame
    """

    endpoint = "https://api.gdc.cancer.gov/files"
    fields = ",".join(fields)
    params = {"filters": filters, "fields": fields, "format": "TSV", "size": n_records}

    print(params)

    resp = requests.post(endpoint, json=params, timeout=1000)
    resp_table = pd.read_table(StringIO(resp.text))

    resp_table = resp_table[resp_table.platform != "Illumina Human Methylation 27"]
    resp_table.columns = [name.split(".")[-1] for name in resp_table.columns]

    return resp_table


def constrain(sample_sheet: pd.DataFrame, n: int) -> pd.DataFrame:
    """

    :param sample_sheet:
    :param n:
    :return: pd.DataFrame
    """
    selected = set()

    for sample in sample_sheet.case_id:
        record = sample_sheet[sample_sheet.case_id == sample]
        record = record["experimental_strategy"].tolist()

        if len(record) == 2 and set(record) == {"Methylation Array", "RNA-Seq"}:
            selected = selected | {sample}

        if len(selected) >= n:
            break

    return sample_sheet[sample_sheet.case_id.isin(selected)]


def build_manifest(samples: list) -> pd.DataFrame:
    """

    :param samples:
    :return: pd.DataFrame
    """
    params = {"ids": samples}
    endpoint = "https://api.gdc.cancer.gov/manifest/"

    resp = requests.post(
        endpoint,
        data=json.dumps(params),
        headers={"Content-Type": "application/json"},
        timeout=1000,
    )

    resp = StringIO(resp.text)
    manifest = pd.read_table(resp)
    return manifest
