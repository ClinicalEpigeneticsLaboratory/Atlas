from os import makedirs
from os.path import join, exists
from subprocess import call
from glob import glob

import pandas as pd

from src.build_manifest_flow import request_gdc_service, build_manifest, constrain
from src.data_proceesing_flow import prepare_exp_data, prepare_meth_data, validate_data, mark_outliers
from src.config import load_config, export_config, add_filters


class Atlas:
    def __init__(
        self,
        project_directory: str,
        primary_site: str,
        diagnosis: str,
        n: int = 200,
        config: str = "config.json",
    ):
        self.project_directory = project_directory
        self.primary_site = primary_site
        self.diagnosis = diagnosis
        self.n = n
        self.config = config

    def prepare_project_dir(self) -> None:
        for sub_dir in ["data/raw", "data/interim", "data/processed", "results"]:
            makedirs(join(self.project_directory, sub_dir), exist_ok=True)

    def prepare_local_config(self):
        if isinstance(self.config, str):
            config_file = load_config(self.config)
            self.config = add_filters(self.primary_site, self.diagnosis, config_file)
            export_config(self.config, join(self.project_directory, "config.json"))

        else:
            print("Config object already created.")

    def create_sample_sheet(self) -> pd.DataFrame:
        raw_sample_sheet = request_gdc_service(
            self.config["FIELDS"], self.config["FILTERS"]
        )

        sample_sheet = constrain(raw_sample_sheet, self.n)
        return sample_sheet

    def create_manifest(self, sample_sheet: pd.DataFrame) -> None:
        manifest_path = join(self.project_directory, "manifest.tsv")

        if exists(manifest_path):
            print("Manifest already exists. Skipping.")
            return None

        manifest = build_manifest(sample_sheet.id.tolist())
        manifest.to_csv(
            manifest_path,
            sep="\t",
            index=False,
        )

    def download_data(
        self,
        manifest_path: str,
        gdc_executable: str,
        n_process: int = 10,
    ) -> None:
        output_path = join(self.project_directory, "data/raw")

        files_in_path = glob(join(output_path, "*"))
        files_declared_in_manifest = pd.read_table(join(self.project_directory, "manifest")).index

        if len(files_in_path) == len(files_declared_in_manifest):
            print("Raw data already exist. Skipping.")
            return None

        command = [
            f"{gdc_executable}",
            "download",
            "-n",
            f"{n_process}",
            "--wait-time",
            "10",
            "-m",
            f"{manifest_path}",
            "-d",
            f"{output_path}",
            "--retry-amount",
            "30",
        ]

        call(command)

    def build_frames(self):
        sample_sheet = pd.read_csv(join(self.project_directory, "sample_sheet.csv"), index_col=0)

        methylation = prepare_meth_data(sample_sheet, join(self.project_directory, "data/raw"))
        expression = prepare_exp_data(sample_sheet, join(self.project_directory, "data/raw"))
        validate_data(methylation, expression)

        met_olist = mark_outliers(methylation.T, join(self.project_directory, "met_pca.png"))
        exp_olist = mark_outliers(expression.T, join(self.project_directory, "met_pca.png"))

        met_to_drop = met_olist[met_olist == "outlier"].index
        exp_to_drop = exp_olist[exp_olist == "outlier"].index

        methylation = methylation.drop(met_to_drop, axis=1)
        expression = expression.drop(exp_to_drop, axis=1)

        methylation.to_parquet(join(self.project_directory, "methylation.parquet"))
        expression.to_parquet(join(self.project_directory, "expression.parquet"))
