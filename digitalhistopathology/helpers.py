#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import argparse
import base64
import io
import json
import os
from pathlib import Path
from typing import Optional, Union

import bokeh.models as bmo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdf2image
import rpy2.robjects as robjects
import scanpy as sc
import scipy
from anndata import AnnData
from bokeh.models import ColumnDataSource
from bokeh.palettes import Category10, Category20
from bokeh.plotting import figure, output_file, save, show
from matplotlib.image import imread
from PIL import Image
from rpy2.robjects import pandas2ri, r
from scipy.cluster.hierarchy import dendrogram
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture


def bubbleHeatmap(
    df_color,
    df_size,
    saving_plot_path="../results/clustering/bubbleheatmap_generatio_pval.pdf",
    height=25,
    width=10,
    plot=True,
    title="",
):
    """Adaptation of the bubble heatmap from R package bubbleHeatmap. df_color and df_size must have the same shape. Be careful the pdf format does not store the
    colorbar, it is a bubbleHeatmap library problem I think...

    Args:
        df_color (pd.DataFrame): DataFrame with values used to defined the color points in the heatmap.
        df_size (_type_): DataFrame with values used to defined the size of the points in the heatmap.
        saving_plot_path (str, optional): Saving path for the plot. Defaults to "../results/clustering/bubbleheatmap_generatio_pval.pdf".
        height (int, optional): Height of the plot. Defaults to 25.
        width (int, optional): Width of the plot. Defaults to 10.
        plot (bool, optional): If we output the plot, for jupyter notebook for example. Defaults to True.
        title (str, optional): Plot title. Defaults to "".
    """
    # TODO: Problem with the pvalue legend color ?? I think problem with pdf saving in R, with png it is ok but the quality is bad.
    df_color = df_color.fillna(999.0)
    df_size = df_size.fillna(999.0)

    pandas2ri.activate()
    df_color_r = pandas2ri.py2rpy(df_color)
    df_size_r = pandas2ri.py2rpy(df_size)

    # Install and load the bubbleHeatmap package
    robjects.r('install.packages("bubbleHeatmap", repos="http://cran.us.r-project.org", quiet=TRUE)')
    robjects.r('library("bubbleHeatmap")')
    # robjects.r('install.packages("colorspace", repos="http://cran.us.r-project.org", quiet=TRUE)')
    # robjects.r('library("colorspace")')

    # Sample R code
    r_code = """
    # Your R code here
    # Make sure to adjust the R code according to your requirements

    # Assuming df1 and df2 are already defined in the R environment
    df_color_r <- as.data.frame({})
    df_size_r <- as.data.frame({})
    df_color_r[df_color_r==999] <- NA
    df_size_r[df_size_r==999] <- NA
    df_color_r <- as.matrix(df_color_r)
    df_size_r <- as.matrix(df_size_r)
    # Create the bubble heatmap
    tree <- bubbleHeatmap(colorMat=df_color_r, sizeMat=df_size_r, treeName="example",
                        leftLabelsTitle=F, showRowBracket=F,
                        rowTitle="Clusters", showColBracket=F, colTitle=F,
                        plotTitle="{}",
                        xTitle=F, yTitle=F,
                        legendTitles=c("Fraction of genes", "-log10(p_values)"))
    """.format(
        df_color_r.r_repr(),
        df_size_r.r_repr(),
        title,
    )

    # Execute the R code
    robjects.r(r_code)

    if saving_plot_path.endswith("pdf"):
        r(f'pdf("{saving_plot_path}", height={height}, width={width})')
    else:
        r(f'png("{saving_plot_path}", height={height*100}, width={width*100})')
    r("grid.draw(tree)")
    r("dev.off()")

    if plot:
        if saving_plot_path.endswith("pdf"):
            image = pdf2image.convert_from_path(saving_plot_path)[0]
        else:
            image = plt.imread(saving_plot_path)
        plt.figure(figsize=(width, height))
        plt.imshow(image)
        plt.axis("off")
        plt.show()


def PIL_image_to_base64(path):
    """
    Take as input the index of an image in our dataset and converts the PIL Image to base64.
    It resturns a url associated to the encoded image.
    """
    im = Image.open(path).convert("RGB")
    buffer = io.BytesIO()
    im.save(buffer, format="jpeg")
    encoded_image = base64.b64encode(buffer.getvalue()).decode()
    im_url = "data:image/jpeg;base64, " + encoded_image
    return im_url


def plot_bokeh(
    df,
    x,
    y,
    img_path="path",
    other_hover_data=["name", "label"],
    title="",
    color_col="name_origin",
    palette=None,
    saving_path=None,
    flip_y_axis=False,
):
    """Bokeh scatterplot with images as hover data.

    Args:
        df (pd.DataFrame): Dataframe with samples as rows and features as columns.
        x (str): Column of df to be plotted in the x axis.
        y (str): Column of df to be plotted in the y axis.
        img_path (str, optional): Column that contains the path to images. Defaults to "path".
        other_hover_data (list, optional): List of columns of df to be included in the hover data. Defaults to ["name", "label"].
        title (str, optional): Title of the plot. Defaults to "".
        color_col (str, optional): _description_. Defaults to "name_origin".
        palette (dict, optional): Palette to color the points, key are the label, values are the colors. Defaults to None.
        saving_path (str, optional): Saving html path. If None, it shows the plot. Defaults to None.
        flip_y_axis (bool, optional): Flip the y axis. Defaults to False.
    """
    df[img_path] = df[img_path].apply(lambda p: PIL_image_to_base64(p)).values
    source_train = ColumnDataSource(df)
    TOOLTIPS = """
    <div>
        <div>
            <img
                src="@{}" height="150" alt="@{}" width="150"
                style="float: left; margin: 0px 15px 15px 0px;"
                border="2"
            ></img>
        </div>
    <div>""".format(
        img_path, img_path
    )
    for hover_data in other_hover_data:
        TOOLTIPS = (
            TOOLTIPS
            + """
            <span style="font-size: 17px; font-weight: bold;">{}</span>
            <span style="font-size: 17px; ">@{}</span>
            </div>
            <div>
        """.format(
                hover_data.capitalize(), hover_data
            )
        )

    fig = figure(
        width=800,
        height=800,
        tooltips=TOOLTIPS,
        title=title,
        tools="pan,wheel_zoom,box_zoom,reset, undo, redo",
        x_axis_label=x,
        y_axis_label=y,
    )

    if color_col is not None:
        factors = list(df[color_col].unique()) if palette is None else list(palette.keys())

        try:
            color_mapper = bmo.CategoricalColorMapper(
                factors=factors, palette=Category10[len(factors)] if palette is None else list(palette.values())
            )
        except Exception as e:
            print(e)
            color_mapper = bmo.CategoricalColorMapper(
                factors=factors, palette=Category20[len(factors)] if palette is None else list(palette.values())
            )

    fig.scatter(
        x,
        y,
        size=8,
        source=source_train,
        alpha=0.6,
        fill_color=None if color_col is None else {"field": color_col, "transform": color_mapper},
        line_alpha=0.3,
    )
    if flip_y_axis:
        fig.y_range.flipped = True

    if saving_path is not None:
        output_file(filename=saving_path, title="Static HTML file")
        save(fig)
    else:
        show(fig)


def intersection(lst1, lst2):
    """Intersection between two lists

    Args:
        lst1 (list): List 1.
        lst2 (list): List 2.

    Returns:
        list: Intersection list
    """
    return sorted(list(set(lst1) & set(list(lst2))))


def corr2_coeff(A, B):
    """From https://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays/30143754#30143754.
    Perform pearson correlation between two matrices efficiently. A and B have the same number of columns.

    Args:
        A (np.array): N x T matrix
        B (np.array): M x T matrix

    Returns:
        np.array: N x M correlation matrix
    """
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None], ssB[None]))


def gmm_threshold(
    data,
    quantile_threshold=0.9,
    background_quantile=True,
    intersection=True,
    density=False,
    bins=50,
    baysian_gmm=True,
    scaling_to_max=True,
    show_plot=True,
    saving_folder=None,
    data_name="",
    title="",
):
    """
    Fit two gaussian distributions to your data and calculated a threshold from the quantile of one of the two gaussians.

    Args:
        data (array): 1D numpy array from which you want to fit two gaussian distributions
        quantile_threshold (float, optional): Quantile threshold you want to apply on the foreground or background gaussian. Defaults to 0.9.
        background_quantile (bool, optional): Apply the quantile threshold on the background (the gaussian with the lowest mean). Defaults to True.
        intersection (bool, optional): To plot or not the intersection between the two gaussians or not. Defaults to True.
        density (bool, optional): To use density instead of count for the histogram. Defaults to False.
        bins (int, optional): Number of bins for the histogram. Defaults to 50.
        baysian_gmm (bool, optional): To use BayesianGaussianMixture instead of GaussianMixture. Defaults to True.
        scaling_to_max (bool, optional): To scale the gaussian to their max inside [mean-std; mean+std], if False it is scaled to the mean inside the same interval. Defaults to True.
        show_plot (bool, optional): To show the plot or not. Defaults to True.
        saving_folder (str, optional): Saving folder to save the figure. Defaults to None.
        data_name (str, optional): Name of your data, to label x axis. Defaults to "".
        title (str, optional): Title of your plot. Defaults to "".

    Returns:
        float: data threshold calculated from the guassian quantile
    """
    if baysian_gmm:
        gmm = BayesianGaussianMixture(n_components=2, random_state=0, max_iter=1000)
    else:
        gmm = GaussianMixture(n_components=2, random_state=0, max_iter=500)

    gmm.fit(np.reshape(data, (data.shape[0], 1)))

    fig, ax = plt.subplots(figsize=(6, 5))

    y_hist, x_hist, _ = plt.hist(data, density=density, alpha=0.5, color="gray", bins=bins)
    x_hist = (x_hist[1:] + x_hist[:-1]) / 2
    hist_df = pd.DataFrame(np.array([y_hist, x_hist]).T, columns=["y", "x"])

    x = np.linspace(data.min(), data.max(), 1000).reshape(-1, 1)

    if gmm.means_[0][0] > gmm.means_[1][0]:
        ind_list = [1, 0]
    else:
        ind_list = [0, 1]

    pdfs = {}
    for i, ind in enumerate(ind_list):
        mean = gmm.means_[ind][0]
        covariance = gmm.covariances_[ind][0, 0]
        dist = scipy.stats.norm(loc=mean, scale=np.sqrt(covariance))
        pdf = dist.pdf(x)

        if (background_quantile and i == 0) or (not background_quantile and i == 1):
            threshold = dist.ppf(quantile_threshold)
            color = "orange"
        else:
            color = "green"

        gauss_max = np.max(pdf)
        hist_max = hist_df[
            np.logical_and(hist_df["x"] > mean - np.sqrt(covariance), hist_df["x"] < mean + np.sqrt(covariance))
        ]["y"]
        if scaling_to_max:
            hist_max = hist_max.max()
        else:
            hist_max = np.nanmean(hist_max.apply(lambda y: float("nan") if y == 0 else y))
        scaling_factor = hist_max / gauss_max
        plt.plot(x, pdf * scaling_factor, label=f"Component {i + 1}", color=color)
        pdfs[f"Component {i + 1}"] = np.reshape(pdf * scaling_factor, -1)

    plt.axvline(threshold, color="red", label=f"{quantile_threshold*100:.0f}th: {threshold:.3f}")
    print("Threshold based on the {}th quantile: {}".format(quantile_threshold * 100, np.round(threshold, 3)))

    if intersection:
        intersection_ind = np.argwhere(np.diff(np.sign(pdfs["Component 1"] - pdfs["Component 2"]))).flatten()
        intersection_point = x[intersection_ind].flatten()
        plt.scatter(
            x[intersection_ind],
            pdfs["Component 1"][intersection_ind],
            c="red",
            label="Intersection: {}".format(np.round(intersection_point.max(), 3)),
            marker="X",
        )
        print("Intersection x point: {}".format(intersection_point[0]))

    plt.legend()
    plt.xlabel(data_name)
    if density:
        plt.ylabel("Density")
    else:
        plt.ylabel("Count")

    plt.title(title)

    if saving_folder is not None:
        plt.savefig(os.path.join(saving_folder, "gmm_threshold_{}.pdf".format(data_name)), bbox_inches="tight")
    if show_plot:
        plt.show()
    else:
        plt.gcf().set_visible(False)

    return threshold


def read_visium(
    path: Union[str, Path],
    genome: Optional[str] = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: Optional[str] = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None,
    image_info_folder: str = "spatial/",
) -> AnnData:
    """\
    From Scanpy library, need to add some changes to be able to use it with my data
    Read 10x-Genomics-formatted visium dataset.

    In addition to reading regular 10x output,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    Based on the `Space Ranger output docs`_.

    See :func:`~scanpy.pl.spatial` for a compatible plotting function.

    .. _Space Ranger output docs: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview

    Parameters
    ----------
    path
        Path to directory for visium datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.

    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:

    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.uns`\\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.
    """
    path = Path(path)
    adata = sc.read_10x_h5(path / count_file, genome=genome)

    adata.uns["spatial"] = dict()

    from h5py import File

    with File(path / count_file, mode="r") as f:
        attrs = dict(f.attrs)
    if library_id is None:
        library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / (image_info_folder + "tissue_positions_list.txt"),
            scalefactors_json_file=path / (image_info_folder + "scalefactors_json.json"),
            hires_image=path / (image_info_folder + "tissue_hires_image.png"),
            lowres_image=path / (image_info_folder + "tissue_lowres_image.png"),
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    logg.warning(f"You seem to be missing an image file.\n" f"Could not find '{f}'.")
                else:
                    raise OSError(f"Could not find '{f}'")

        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["hires", "lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(str(files[f"{res}_image"]))
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(files["scalefactors_json_file"].read_bytes())

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        # read coordinates
        positions = pd.read_csv(files["tissue_positions_file"], header=None)
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        positions.index = positions["barcode"]

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
        adata.obs.drop(
            columns=["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(source_image_path)

    return adata


def parsing_oldST(
    adata: AnnData,
    coordinates_file: Union[Path, str],
    copy: bool = True,
) -> Optional[AnnData]:
    """\
    From stlearn/adds/parsing.py, need to comment the line in which they filter genes !
    https://github.com/BiomedicalMachineLearning/stLearn/blob/master/stlearn/adds/parsing.py
    Parsing the old spaital transcriptomics data.

    Parameters
    ----------
    adata
        Annotated data matrix.
    coordinates_file
        Coordinate file generated by st_spot_detector.
    copy
        Return a copy instead of writing to adata.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    **imagecol** and **imagerow** : `adata.obs` field
        Spatial information of the tissue image.
    """

    # Get a map of the new coordinates
    new_coordinates = dict()
    with open(coordinates_file, "r") as filehandler:
        for line in filehandler.readlines():
            tokens = line.split(sep=",")[1:] if "," in line else line.split(sep="\t")
            assert len(tokens) >= 6 or len(tokens) == 4
            if tokens[0] != "x":
                old_x = int(tokens[0])
                old_y = int(tokens[1])
                new_x = round(float(tokens[2]), 2)
                new_y = round(float(tokens[3]), 2)
                if len(tokens) >= 6:
                    pixel_x = float(tokens[4])
                    pixel_y = float(tokens[5])
                    new_coordinates[(old_x, old_y)] = (pixel_x, pixel_y)
                else:
                    raise ValueError(
                        "Error, output format is pixel coordinates but\n " "the coordinates file only contains 4 columns\n"
                    )

    counts_table = adata.to_df()
    new_index_values = list()

    imgcol = []
    imgrow = []
    for index in counts_table.index:
        tokens = index.split("x")
        x = int(tokens[0])
        y = int(tokens[1])
        try:
            new_x, new_y = new_coordinates[(x, y)]
            imgcol.append(new_x)
            imgrow.append(new_y)

            new_index_values.append("{0}x{1}".format(new_x, new_y))
        except KeyError:
            counts_table.drop(index, inplace=True)

    # Assign the new indexes
    # counts_table.index = new_index_values

    # Remove genes that have now a total count of zero
    # counts_table = counts_table.transpose()[counts_table.sum(axis=0) > 0].transpose()

    adata = AnnData(counts_table)

    adata.obs["imagecol"] = imgcol
    adata.obs["imagerow"] = imgrow

    adata.obsm["spatial"] = np.c_[[imgcol, imgrow]].reshape(-1, 2)

    return adata if copy else None


def str2bool(v):
    """Change str to bool.

    Args:
        v (str): Parse argument.

    Raises:
        argparse.ArgumentTypeError

    Returns:
        bool: True or False
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def quantile_normalize(df: pd.DataFrame):
    """Performs quantile normalization across the columns. This function uses the function normalizeQuantiles from R.

    Args:
        df: pd.DataFrame on which quantile normalization will be performed on the columns.

    Returns:
        pd.DataFrame: The normalized dataframe.
    """

    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr

    with (robjects.default_converter + pandas2ri.converter).context():
        r_bulk = robjects.conversion.get_conversion().py2rpy(df)

    limma = importr("limma")
    r_bulk_norm = limma.normalizeQuantiles(r_bulk)
    df_norm = pd.DataFrame(dict(zip(r_bulk_norm.names, list(r_bulk_norm))))
    df_norm.index = df.index
    return df_norm




class NumpyEncoder(json.JSONEncoder):
    """Transform array to list while saving json. From: https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable."""

    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
