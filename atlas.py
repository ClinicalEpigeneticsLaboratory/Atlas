from os import makedirs
from os.path import join, exists
from subprocess import call

import pandas as pd

from src.build_manifest_flow import request_gdc_service, build_manifest, constrain
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
        command = [
            f"{gdc_executable}",
            "download",
            "-n",
            f"{n_process}",
            "--wait-time",
            "5",
            "-m",
            f"{manifest_path}",
            "-d",
            f"{output_path}",
            "--retry-amount",
            "25",
        ]

        call(command)
