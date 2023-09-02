from tqdm import tqdm
from glob import glob
from os.path import join
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.neighbors import LocalOutlierFactor
from sklearn.preprocessing import StandardScaler


class LowFractionOfCommonSamplesError(Exception):
    pass


def prepare_meth_data(sample_sheet: pd.DataFrame, data_source: str, na_limit: float = 0.2):
    dataset = []
    for file in tqdm(sample_sheet.id):
        path = join(data_source, file, "*level3betas.txt")
        path = glob(path)

        if not len(path):
            continue

        path = path[0]
        frame = pd.read_table(path, index_col=0, header=None)
        frame.columns = [sample_sheet[sample_sheet.id == file]["case_id"].squeeze()]
        frame.index.name = ""

        if frame.iloc[:, 0].isna().sum() / len(frame) > na_limit:
            print(f"Skipping file {file} due to large NaN frequency > {na_limit}")
            continue

        dataset.append(frame)

    dataset = pd.concat(dataset, axis=1)
    dataset = dataset.dropna()
    dataset = dataset.loc[[name for name in dataset.index if str(name).startswith("cg")]]
    return dataset


def prepare_exp_data(sample_sheet: pd.DataFrame, data_source: str, median_tpm_limit: int = 1):
    dataset = []
    for file in tqdm(sample_sheet.id):
        path = join(data_source, file, "*star_gene_counts.tsv")
        path = glob(path)

        if not len(path):
            continue

        path = path[0]
        frame = pd.read_table(path, index_col=0, comment="#")

        frame = frame[["gene_name", "tpm_unstranded"]]
        frame = frame.set_index("gene_name")

        frame.columns = [sample_sheet[sample_sheet.id == file]["case_id"].squeeze()]
        frame.index.name = ""
        dataset.append(frame)

    dataset = pd.concat(dataset, axis=1)
    dataset = dataset.dropna()

    median = dataset.median(axis=1)
    to_drop = median[median < median_tpm_limit].index
    dataset = dataset.drop(to_drop)

    return dataset


def validate_data(meth_dataset: pd.DataFrame, exp_dataset: pd.DataFrame) -> None:
    intersect = set(meth_dataset.columns).intersection(set(exp_dataset.columns))

    if len(intersect) < 100:
        raise LowFractionOfCommonSamplesError


def mark_outliers(raw_dataset: pd.DataFrame, pca_plot_path: str, n: int = 25) -> pd.Series:
    pca = PCA(n_components=n)

    dataset = StandardScaler().fit_transform(raw_dataset)
    dataset = pca.fit_transform(dataset)
    dataset = pd.DataFrame(dataset, index=raw_dataset.index, columns=[f"PC {i + 1}" for i in range(n)])

    lof = LocalOutlierFactor()
    prediction = lof.fit_predict(dataset)

    outliers = pd.Series(prediction, index=raw_dataset.index, name="Outliers")
    outliers = outliers.replace({1: "inlier", -1: "outlier"})

    pc1, pc2 = pca.explained_variance_ratio_[:2]
    pc1, pc2 = round(pc1 * 100, 1), round(pc2 * 100, 1)
    pc1_label = f"PC 1 {pc1}%"
    pc2_label = f"PC 2 {pc2}%"

    fig = px.scatter(dataset, x="PC 1", y="PC 2", color=outliers, labels={"PC 1": pc1_label, "PC 2": pc2_label})
    fig.update_layout(width=600, height=600, font={"size": 16}, legend={"title": ""})
    fig.write_image(pca_plot_path)

    return outliers
