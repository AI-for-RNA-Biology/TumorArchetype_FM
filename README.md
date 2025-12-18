# Image-Based Discovery of a Recurrent Tumor Archetype Characterized by Aberrant RNA Splicing and Associated with Poor Survival in Breast Cancer

## Overview
Here, we develop a pipeline to systematically evaluate the biological concepts encoded within histopathology foundation models (hFMs) using molecular data. We also perform extended pre-training of [UNI](https://github.com/mahmoodlab/UNI) to identify optimal conditions that enhance the model’s ability to encode richer, tumor tissue-specific biological concepts.

<img width="4961" height="7016" alt="Figure 7" src="https://github.com/user-attachments/assets/517f1f6c-9a9a-4272-91dd-919ab01a36b5" />

This code was developed on a 64-bit Debian GNU/Linux system (x86_64), running kernel version 6.1.0-26-amd64 (Debian 6.1.112-1), with dynamic preemption enabled (PREEMPT_DYNAMIC). The code to perform extended pretraining is available in the [dinov2](https://github.com/AI-for-RNA-Biology/dinov2/tree/19a441962b97527abc3a8ca2da8a2938aa7c663e) folder, while the code to assess compute and evaluate the resulting embeddings is located in the [digitalhistopathology](./digitalhistopathology/) and [scripts](./scripts/) folders.

## 1. Installation

Clone the repository, and initialize the submodules:
```bash
git clone git@github.com:AI-for-RNA-Biology/TumorArchetype_FM.git # using SSH
cd TumorArchetype_FM
git submodule update --init --recursive
```

The environment can be installed via conda:
```bash
conda env create -f environment.yaml
```

After the installation, you can activate the environment with:
```bash
conda activate digitalhisto
```

Then add the following libraries with the commands:
```bash
pip install stlearn --no-deps --force-reinstall
pip install graphviz --no-deps --force-reinstall
```

Install R packages:
```bash
R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("EBImage")'
```


**BE CAREFUL: If you want to run CTransPath you will need a specific version of the timm library.**

Go to their [Github page](https://github.com/Xiyue-Wang/TransPath) under 1.CTransPath and install the modified timm library on top of the anaconda environment.
When you want to run other foundation models like Uni or ProvGigaPath, override the timm library with version 1.0.7:

```bash
pip install timm==1.0.7
```

If you want to be safer, create two anaconda environements starting from the yaml file with one specific to CTransPath and one for the other models.

The typical installation time in around 25 minutes.

## 2. Data

### HER2-positive breast cancer dataset [[1](#ref1)]

The data can be downloaded from [Zenodo](https://zenodo.org/records/4751624).

- Rename the principal folder as `HER2`
- Unzip all the folders: `count-matrices.zip`, `images.zip`, `meta.zip` and `spot-selections.zip`
- Place the principal folder in `data/`

### Triple negative breast cancer dataset [[2](#ref10)]

The data can be downloaded from [Zenodo](https://zenodo.org/records/8135722).

- Rename the principal folder as `TNBC`
- Unzip all the folders
- Place the principal folder in `data/`

All other new datasets can be added by creating a new class in `digitalhistopathology/datasets/real_datasets.py`

## 3. HuggingFace Hub

UNI and Prov-GigaPath are pretrained models whose access must be granted before using them. You must have a hugging face account, agree to the outlined terms of use, and request access.

Your hugging face access tokens (must remain private) should be set in a new file `digitalhistopathology/access_token.py` with the following content:
```python
READ_ONLY_HF_TOKEN = "hf_xxxxxxxxxxxxxxxxxxx"
```
The first time loading the model, the weights will be downloaded to the hugging face hub cache in your home directory (`~/.cache/huggingface/hub`).

More information on [cache management](https://huggingface.co/docs/datasets/en/cache) and on [environment variables](https://huggingface.co/docs/huggingface_hub/en/package_reference/environment_variables).

## 4. Pretrained models

For SimCLR and CTransPath place the pretrained model weights you want to test in the `pretrained_models` folder. All other new pretrained models can be added by creating a new class in `digitalhistopathology/models.py`

### 4.1 SimCLR [[3](#ref3)]

The weigths can be found to the [github page](https://github.com/ozanciga/self-supervised-histopathology) of the projet.

### 4.2 CTransPath [[4](#ref4)]

The weigths can be found to the [github page](https://github.com/Xiyue-Wang/TransPath) of the projet.

### 4.3 UNI [[5](#ref5)]

The weights can be found on the [huggingface page](https://huggingface.co/MahmoodLab/UNI) of the project.

### 4.4 Prov-GigaPath [[6](#ref6)]

The weights can be found on the [huggingface page](https://huggingface.co/prov-gigapath/prov-gigapath) of the project.

### 4.5 Virchow2 [[7](#ref7)]

The weights can be found on the [huggingface page](https://huggingface.co/paige-ai/Virchow2) of the project.

### 4.6 UNI2 [[8](#ref)]

The weights can be found on the [huggingface page](https://huggingface.co/MahmoodLab/UNI2-h) of the project.

## 5. Full pipeline

The following pipeline aims at reproducing the results from the original paper. 

### Description

Here are the different steps of the pipeline that needs a dataset, a model and a name to be initialized.

### 5.1. Compute patches:

First, we need to compute patches from the desired dataset. This can be done using the script `scripts/compute_patches.py`, which takes one argument:
- `--dataset`: Dataset to process. Choose among "HER2" and "TNBC". Default is "HER2".

Example usage: 
```bash
python scripts/compute_patches.py --dataset "HER2"
```
The resulting patches will be produced in the `results/{dataset}/compute_patches/all/` folder. 

__Note: For HER2 dataset, patient A has been excluded from the analysis for suspicion of technical bias.__

### 5.2. Compute image embeddings and k-NN selection of invasive cancer patches:

This step is made to compute the embeddings using a pre-trained model, and to apply k-NN on them to extract the patches classified as "invasive cancer patches". The shannon entropy computed on the explained variance distribution of SVD components is also computed, and molecular data are formatted for later use. This is done using the script `scripts/pipeline.py`. 

Here are the arguments you need to provide:
- `--model_name`: Pretrained model to be used. The supported values are "uni", "ctranspath", "provgigapath", "simclr", "virchow". You can add other models in `digitalhistopathology/models.py`. Default is "uni".
- `--retrained_model_path`: If you retrained the pre-trained model and want to extract the associated embeddings, specifiy the path to the weights file. 
- `--patches_folder`: Folder to the patches dataset you want to compute the embeddings from. Default is `results/HER2/compute_patches/all/`.
- `--pipeline_name`: How you want to name the pipeline
- `--results_folder`: Folder to save the embeddings. Default is `results/`.
- `--dataset`: Dataset name. Choose among "HER2" and "TNBC". Default is "HER2".
All results will be saved at `{results_folder}/{dataset}/pipeline/{pipeline_name}/`.

In more details here are the steps:
- Patches image embeddings will be computed and saved as `{results_folder}/{dataset}/pipeline/{pipeline_name}/image_embedding.h5ad`:
- Select invasive cancer patches by computing KNN on `image_embeddding` by using label column as training data. At the end, we selected only invasive_cancer label to go more in depth into tumor heterogeneity.
   - Boxplot of the F1 and accuracy across patient to choose the optimal k of knn and barplots with patient and label fraction in each knn cluster are saved in the folder `{results_folder}/{dataset}/pipeline/{pipeline_name}/select_invasive_cancer`.
   - The invasive images embedding is saved under `{results_folder}/{dataset}/pipeline/{pipeline_name}/invasive_image_embedding.h5ad`.
- Shannon entropy is computed on the SVD decomposition and results are saved in the folder `{results_folder}/{dataset}/pipeline/{pipeline_name}/shannon_entropy`.

Example usage to extract embeddings from UNI using HER2 dataset:
```bash
python scripts/pipeline.py --model_name uni --patches_folder results/HER2/compute_patches/all/ --pipeline_name uni --results_folder results --dataset HER2
```

Note: If you just want to extract embeddings without computing the shannon entropy or running the k-NN re-annotation, you can use the ready-to-use script `scripts/compute_embeddings.py` that takes the same arguments. 

This step requires a GPU. It was run on a rtx3090 GPU, 12 CPUs and 128GB of RAM. You can also compute without GPU as it is useful only for inference. It lasts around 20 minutes.

The config files for this step and all models are available under `config/pipeline/base_models`. You can refer to them to know the parameters you need to pass to the function.

### 5.3. Create the folder for invasive cancer patches only:

This step enables to create a new patches dataset, subset of the original dataset and filtered using a csv file. It is done with the script `scripts/create_subset_hdf5.py` that takes three arguments:

- `--original_hdf5`: Path to the original HDF5 file.
- `--csv_path`: Path to the CSV file containing the names of the patches to include.
- `--new_hdf5`: Path to save the new HDF5 file.


In the context of our study, we were interested in the invasive cancer patches, and after the kNN-relabeling of step 2, we created a filtered dataset to keep only the patches relabeled as "invasive cancer" from the pipeline applied to the embeddings.

Example usage with UNI model and HER2 dataset:
```bash
python scripts/create_subset_hdf5.py \
--original_hdf5 results/HER2/compute_patches/all/patches.hdf5 \
--csv_path results/HER2/pipeline/uni/select_invasive_cancer/selected_invasive_cancer.csv \
--new_hdf5 results/HER2/compute_patches/invasive/patches.hdf5
``` 

This step is really fast and can be run on a single CPU with 16GB of memory. 

### 5.4. Extended-pretraining:

**5.4.1. Perform extended-pretraining**:

In this step, we performed extended pre-training of UNI, using the invasive cancer patches selected above. The final models are called T-UNI models. Different T-UNI models have been generated, varying the extended pretraining strategy (full or ExPLoRa), the loss (KDE or KoLeo), and finally, the number of prototypes. All used config files are available in `config/extended_pretraining`. The script to perform extended pre-training is `dinov2/dinov2/train/train.py` and can be used as follows:

```bash
PYTHONPATH=$(pwd)/dinov2 python dinov2/dinov2/train/train.py config/extended_pretraining/{dataset}/{config_file}.yaml
```

For more information about the extended pretraining, please see the [dinov2 README](https://github.com/AI-for-RNA-Biology/dinov2/tree/19a441962b97527abc3a8ca2da8a2938aa7c663e/README.md). Note that you will need another environment to run it. The [model card](./dinov2/MODEL_CARD_T_UNI.md) is available for the T-UNI models. You can also download direclty the model weights on [Zenodo](https://zenodo.org/records/15053890) and place them in the folder `pretrained_models/{dataset}`. For the rest of the code to run, please put each model in a directory named with the model name, then rename the weights file as `epoch_1.pth`.

Example: 
When you download the model `HER2_uni_full_kde_16384_prototypes.pth` on Zenodo, place it into `pretrained_models/HER2/HER2_uni_full_kde_16384_prototypes/epoch_1.pth`.

**5.4.2 Run the embedding extraction on the extended pre-trained models**:

Example usage to extract embeddings from a retrained version of UNI located in `pretrained_models/HER2/HER2_uni_full_kde_16384_prototypes` under `epoch_1.pth`:

```bash
python scripts/pipeline.py --model_name uni --retrained_model_path pretrained_models/HER2/HER2_uni_full_kde_16384_prototypes/epoch_1.pth --patches_folder results/HER2/compute_patches/ --pipeline_name HER2_uni_full_kde_16384_prototypes --results_folder results/HER2/pipeline 
```

The config files to extract embeddings from all T-UNI models are available under `config/pipeline/T-UNI`.

### 5.5 Handcrafted features extraction

In this step, segmentation of the nuclei has been performed using CellViT [[9](#ref)]. Features related to the morphology and texture of the nuclei have been extracted using scMTOP and correspond to the “Nuclei-Morph”, and “Nuclei-Texture” features respectively. We have enhanced the “Nuclei-Morph” category by computing Zernike moments of each cell using the Mahotas python package. We have adapted scMTOP to add color descriptors of the nuclei, referred to as “Nuclei-Color” features, namely mean, skewness, kurtosis and entropy of each RGB color channel, as well as intensity and transparency. Finally, statistics on the cell type composition have been computed based on cell labels outputted by CellViT, and are referred to as “Nuclei-Composition”. Features related to the extracellular matrix color (“ExtraCell-Color” features) and texture (“ExtraCell-Texture”) have been computed using scikit-image. Similarly, texture and color features at the entire patch level have been extracted and are referred to as “WholePatch-Texture” and “WholePatch-Color”.

**5.5.1 Nuclei segmentation using CellViT**:

The [CellViT](./CellViT/) folder is a Git subtree cloned from the [original repository](https://github.com/TIO-IKIM/CellViT), with only minimal modifications made to adapt it to our environment. 

First, you need to set up the CellViT environment:
```bash
conda env create -f cellvit_env.yml
conda activate cellvit_env
```

Then, download the [CellViT checkpoint] (https://drive.google.com/uc?export=download&id=1wP4WhHLNwyJv97AK42pWK8kPoWlrqi30) and place it in the `CellViT/models/pretrained/` folder. This model works with slide magnifications 20x (used for both HER2 and TNBC), model for 40x is also available [here](https://drive.google.com/uc?export=download&id=1MvRKNzDW2eHbQb5rAgTEp6s2zAXHixRV).


Finally, run the segmentation using the script `digitalhistopathology/engineered_features/cell_segmentor.py`, with the following arguments: 
- `--segmentation_mode`: Segmentation mode to use. For CellViT, set this to "cellvit".
- `--magnification`: Magnification level of the images (e.g., 20x).
- `--mpp`: Microns per pixel for the images.
- `--patches_info_filename`: File containing metadata about image patches.
- `--list_wsi_filenames`: List of whole slide image filenames to process.
- `--dataset_name`: Name of the dataset being processed.
- `--model_path`: Path to the pretrained CellViT model weights.
- `--result_saving_folder`: Directory where the segmentation results will be saved.

_Note: If `--patches_info_filename` is provided, the `--list_wsi_filenames` argument will be ignored. Use `--list_wsi_filenames` when you need to process specific individual files._

Example usage for HER2 dataset:
```bash
PYTHONPATH=$(pwd)/CellViT python3 digitalhistopathology/engineered_features/cell_segmentor.py \
--segmentation_mode "cellvit" \
--magnification 20 \
--mpp 1 \
--patches_info_filename "results/HER2/compute_patches/all/patches_info.pkl.gz" \
--dataset_name "HER2" \
--model_path "CellViT/models/pretrained/" \
--result_saving_folder "results/HER2/"
```

The results will be saved in `results/{dataset}/segmentation/`. This step took around 30 minutes for HER2 dataset on a h100 GPU with 12 cpus per task and 64GB of memory.

**5.5.2 Handcrafted features computation**: 

This step computes handcrafted features such as morphology, texture, and color descriptors of nuclei and extracellular matrix. It also includes Zernike moments for nuclei and cell type composition statistics.

Here are the parameters for the `digitalhistopathology/engineered_features/engineered_features.py` script:
- `--method`: Feature extraction method to use (e.g., morphology, texture, color).
- `--path_to_cellvit_folder`: Path to the folder containing CellViT segmentation results.
- `--result_saving_folder`: Directory where the computed features will be saved.
- `--dataset_name`: Name of the dataset being processed.
- `--patches_info_filename`: File containing metadata about image patches.
- `--list_wsi_filenames`: List of whole slide image filenames to process.
- `--save_individual_wsi`: Flag to save results for each WSI individually.
- `--zernike_per_nuclei`: Flag to compute Zernike moments for each nucleus.
- `--num_cores`: Number of CPU cores to use for parallel processing.

_Note: If `--patches_info_filename` is provided, the `--list_wsi_filenames` argument will be ignored. Use `--list_wsi_filenames` when you need to process specific individual files._

Example for HER2 dataset without Zernike moments:
```bash
PYTHONPATH=$(pwd)/CellViT python3 digitalhistopathology/engineered_features/engineered_features.py \
--method "scMTOP" \
--path_to_cellvit_folder "results/HER2/segmentation/" \
--result_saving_folder "results/HER2/engineered_features/scMTOP/" \
--dataset_name "HER2" \
--patches_info_filename "results/HER2/compute_patches/all/patches_info.pkl.gz" \
--num_cores 16
```

The results will be saved in `results/{dataset}/engineered_features/scMTOP`. This step is particularly long to run, especially due to zernike moments computation. It can indeed take up to 7 hours per slide, with 16 CPUs and 24GB of memory. Therefore we recommand lauching the script in parallel for each WSI. In future versions, this script will include built-in parallelism to streamline the process, and the operations will be optimized in order to save computation time. 


### 5.6 Molecular data preparation

This section describes how to preprocess, load, and analyze molecular data for the HER2 and TNBC datasets. The scripts are available in the [scripts/HER2](scripts/HER2) and [scripts/TNBC](scripts/TNBC) folders.


**5.6.0 Preprocessing (TNBC only)**

This step converts raw RDS count objects into a filtered, patient-aware Seurat object. Skip this section for HER2, which is already preprocessed.

1. Install the required dependency, and set the library path:
    ```bash
    pip install gtfparse=2.5.0
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
    ```

2. Download the [GENCODE v38 GTF file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz), extract it and place it in `annotation/gencode.v38.annotation.gtf`.

3. Run the script `scripts/TNBC/0_1_process_rds_memory_efficient.py` to process the RDS files into a concatenated count matrix:
    ```bash
    python scripts/TNBC/0_1_process_rds_memory_efficient.py \
    --rds-dir "data/TNBC/Robjects/countsNonCorrected" \
    --gtf-path "annotation/gencode.v38.annotation.gtf" \
    --output-dir "results/TNBC/molecular/" \
    --batch-size 10 \
    --output-name "TNBC_all_samples_coding_genes" \
    --save-csv \
    --keep-intermediates
    ```

    This step loads RDS count objects sample-by-sample, keeps protein-coding genes only and writes intermediate parquet batches.


4. Run the script `scripts/TNBC/0_2_process_tnbc_seurat.R` to process patient-wise filtering and create Seurat object:
    ```bash
    Rscript scripts/TNBC/0_2_process_tnbc_seurat.R \
        "parquet_batch" \
        "" \
        "results/TNBC/molecular/intermediate" \
        "results/TNBC/molecular/filtering" \
        0.10 \
        FALSE \
        results/compute_patches/TNBC/spots_labels.csv
    ```

    This step applies patient-wise bimodal filtering, selects genes present in ≥ X% (default: 10% of all genes) of patients, and creates a filtered Seurat object.


**5.6.1 Load molecular data (HER2 and TNBC)**
To load the molecular data properly, run the dataset-specific script `scripts/{dataset}/1_load_molecular.py`. 

Example usage for HER2 dataset:
```bash
python scripts/HER2/1_load_molecular_data.py \
--gene_embedding_saving_path results/HER2/molecular \
--patches_folder results/HER2/compute_patches/all
```

This step is really fast (a few minutes) and can be run on a single cpu wih 16GB.


**5.6.2 Tissue-level differential gene expression (DGE)**
This step identifies biological pathways associated with each tissue type by performing differential gene expression (DGE) analysis on the spatial transcriptomics data.

Example usage for HER2 dataset and UNI model:
```bash 
Rscript scripts/HER2/2_tissue_DGE.R uni
```

Example usage for TNBC dataset:
```bash
Rscript scripts/TNBC/2_tissue_DGE.R \
    "results/TNBC/molecular/filtering/TNBC_filtered_seurat_object.rds" \
    "results/TNBC/molecular/" \
    TRUE
```
The parameter `TRUE` indicates that the genes are upregulated, while `FALSE` indicates that the genes are downregulated.


**5.6.3 DGE between invasive clusters**

This step performs differential analysis between invasive cancer clusters inferred by trained models.

Example usage for HER2 dataset and UNI model:
```bash
Rscript scripts/HER2/3_DGE_between_invasive_clusters.R uni
```
Same command structure for TNBC dataset.

This step loads invasive cancer cluster assignments, computes molecular profiles for each cluster, and performs pathway enrichment analysis.


### 5.7. **Run the benchmarks**:

All the benchmarks can be run usin
g the script `script/benchmark_task.py`. Independent on the benchmark task, you will need to provide the following arguments: 
- `--path_to_pipeline`: The path to the folder(s) containing the pipeline results for the models you want to benchmark. This is the path to the `{results_folder}/{pipeline_name}/` as described in **5.2**. Multiple paths can be provided to process multiple models, separated by spaces.
- `--pipelines_list`: A list of names corresponding to the pipelines provided in `--path_to_pipeline`. These names will be used to label the results.
- `--saving_folder`: The folder where the benchmark results will be saved.
- `--dataset`: The name of the dataset being used for benchmarking. This is used for labeling and organizing results.
- `--engineered_features_saving_folder`: The folder where engineered/hancrafted features, such as those computed from image embeddings, are saved.
- `--extension`: The file format for saving plots and figures generated during the benchmarking process (e.g., `pdf`, `png`).


**5.7.1 Shannon entropy**:
This step computes the shannon entropy of the explained variance distribution from SVD components computed on the raw hFM embeddings (`image_embedding`). It is done using the script `python3 benchmark_task.py --benchmark_task shannon_entropy`.

On top of the arguments described above, used in all benchmark tasks, the shannon entropy task takes:
- `--min_comp`: The number of components on which the shannon entropy will be computed. It should correspond to the smallest embedding size. In our case it was 512, the size of SimCLR embeddings.
- `--pct_variance`: This parameter determines the number of components to keep based on the cumulative variance explained, ensuring that the specified percentage of the total variance is captured.

You can also load the parameters from a config file.

Example usage:

```bash
python scripts/benchmark_task.py --config config/benchmarks/benchmark_shannon_base_models.json
```

For a list of 7 models, this step should take around 30 minutes on a single cpu using 32GB of memory.

The config files we used in this paper are `config/benchmarks/benchmark_shannon_base_models.json`, `config/benchmarks/benchmark_shannon_uni_explora.json` and `config/benchmarks/benchmark_shannon_uni_full.json`.


**5.7.2 Clustering pipeline with performance assessment using ARI score**:

In this step, clustering was performed on the image representations obtained from various hFMs or from our selected handcrafted features taken as a whole. Performance was measured using  Adjusted Rand Index (ARI) scores computed using the original tissue types as true labels (“ARI with tissue type”). Before applying the k-means clustering algorithm with k set as the known number of classes, 2-dimensional UMAP embeddings were computed to account for non-linear relations between image representations, validating various hyperparameters: n_neighbors ranging from 10 to 400, and min_dist values of 0.001 and 0.1. For each configuration, hyperparameters maximizing the  ARI with tissue type score were selected. Clustering is performed either on the whole dataset or for each patient individually. The UMAP-k-means approach was compared to raw-k-means where k-means is directly computed on the raw image representations, and to SVD5-kmeans which factor image representations using  the first 5 principal components of SVD before performing k-means. This step is done using the script `python3 benchmark_task.py --benchmark_task unsupervised_clustering_ARI`.

On top of the arguments described above, you will need to provide a clustering algorithm with the argument `--clustering_algo`. In our analysis, it was set to `kmeans`.

For a single model, this step takes approximately 3 hours and 30 minutes on a single CPU with 24 GB of memory. The computation time is primarily due to the numerous UMAP projections performed across different parameter combinations. We therefore recommend running it in parallel for each model (e.g., launch separated jobs) to generate the necessary files efficiently. Once all individual computations are complete, the script can be rerun with the full list of models to generate comparative plots. In future versions, this script will include built-in parallelism to streamline the process.

Example usage:

```bash
python scripts/benchmark_task.py --config config/benchmarks/benchmark_clustering_ARI_base_models.json

```
The config files we used in this paper are `config/benchmarks/benchmark_clustering_ARI_base_models.json`, `config/benchmarks/benchmark_clustering_ARI_uni_explora.json` and `config/benchmarks/benchmark_clustering_ARI_uni_full.json`.

**5.7.3 Linear regression to predict handcrafted features from image embeddings**:

In this step, we predict each handcrafted feature individually from hFM embeddings using linear regression. A 5-fold cross-validation approach is employed, and the mean R² score across the five folds is reported for each feature.

On top of the classic arguments needed for a benchmark task, you will need to provide regression-specific parameters:
- `--regression_type`: Type of regression to use (e.g., `linear` for linear regression).
- `--n_splits`: Number of splits for cross-validation (default is 5).
- `--on_invasive`: Flag to indicate whether to perform regression only on invasive cancer patches.

Example usage:
```bash
python scripts/benchmark_task.py --config config/benchmarks/benchmark_regression_base_models.json
```
And to perform regression on invasive cancer patches only:
```bash
python scripts/benchmark_task.py --config config/benchmarks/benchmark_regression_base_models.json --on_invasive
```

The results will include the R² scores for each handcrafted feature, saved in  `{saving_folder}/regression/{regression_type}`. These scores provide insights into how well the hFM embeddings can predict specific handcrafted features.

The config files we used in this paper are `config/benchmarks/benchmark_regression_base_models.json`, `config/benchmarks/benchmark_regression_uni_explora.json` and `config/benchmarks/benchmark_regression_uni_full.json`.

This step takes approximately 5 hours for a single model on a single cpu with 16GB of memory. This is why we recommand lauching different models in parallel. In future versions, this script will include built-in parallelism to streamline the process.


**5.7.4 Identification of invasive cancer archetypes**:

Clustering was performed with the UMAP-k-means strategy on the invasive cancer patches only. Hyperparameter validation for both UMAP and k-means was performed toward maximizing jointly 1) the silhouette score, and 2) the “batch effect mitigation” score, defined as 1-”ARI with patient”, computed using patient labels instead of tissue labels on patches from invasive cancer. 

The wasserstein distance is then calculated in the image embedding space between the identified clusters. The wasserstein distance is also computed patient-wise between clusters in the molecular space. 
Remark that different hFM representation spaces can have different volumes but similar discernibility between archetypes or clusters. Hence for comparing such spaces in an unbiased manner,  the quantized Wasserstein distance in the image space was computed on normalized representations. Specifically, we computed the maximal diameter α of the UNI representation space by finding the maximal distance between the two most distant points, and normalized the other embedding spaces by multiplying them by their maximal diameters and dividing by α.

On top of the classic arguments needed for a benchmark task, you will need to provide the specific arguments:
- `--ref_model_emb`: path to the reference image embedding, used to compute the normalized quantized wasserstein embeddings. In the paper, we use UNI as reference embedding. 
- `--molecular_emb`: path to the .h5ad molecular embedding you want to use to compute the patient-wise inter-clusters quantized wassestein distance.
- `--molecular_name`: name of the molecular embeddings for labels.
- `--clustering_algo`: clustering algorithm to use for clustering of invasive cancer patches

Example usage:

```bash
python scripts/benchmark_task.py --config config/benchmarks/benchmark_invasive_full_models.json

```

This step can take up to 8 hours on a single cpu with 32GB of memory. Therefore, we recommand running each model individually (i.e, lauching different jobs). The config files we used in this paper are `config/benchmarks/benchmark_invasive_base_models.json`, `config/benchmarks/benchmark_invasive_uni_explora.json` and `config/benchmarks/benchmark_invasive_uni_full.json`.

## 6 Notebooks:

These notebooks analyze the results produced by the full pipeline through additional analyses and vizualizations. You can run them in their order of appearance. All notebooks load the required files from a config file, specific to the dataset: `config/config_notebooks_{dataset}.json`.

**Description of the notebooks:**

Common to both datasets:
- `1_Unsupervised_clustering_base.ipynb`: Analyze the results of the clustering and the ARI scores for the base models.
- `2_Regression_handcrafted_features.ipynb`: Analyze the results from linear regression, for both the base models and the retrained models. 
- `3_Unsupervised_clustering.ipynb`: Analyze the results of the clustering and the ARI scores for the retrained models. 
- `4_Shannon_entropy.ipynb`: Computes the Shannon entropy of the retrained models.
- `5_Invasive_cancer_clustering.ipynb`: Focuses on finding the best UMAP and k-means parameters for clustering invasive cancer patches to identify distinct archetypes. 
- `6_Invasive_clusters_image_wasserstein.ipynb`: Visualizes Wasserstein distances between invasive cancer clusters in the image embedding space. 
- `7_Invasive_clusters_molecular_wassertein.ipynb`: Visualizes Wasserstein distances between invasive cancer clusters in the molecular embedding space. 
- `8_Invasive_clusters_images.ipynb`: Visualizes images from invasive cancer clusters. 
- `9_1_Overall_batch.ipynb`: Analyzes batch effects across all tissue types.
- `9_2_Invasive_batch.ipynb`: Analyzes batch effects in invasive cancer clusters.

Specific to HER2:
- `10_Invasive_clusters_viz.ipynb`: Creates visualizations for invasive cancer clusters, on the original WSIs. 
- `11_Gene_expression_viz.ipynb`: Visualizes gene expression patterns across clusters on the WSIs. 
- `S1_Patient_A.ipynb`: Supplementary analysis showing why patient A was excluded from all analyses.
- `S2_Macenko.ipynb`: Supplementary analysis showing the results obtained using Macenko stain normalization.

Specific to TNBC:
- `10_Invasive_clusters_among_patients.ipynb`: Visualizes invasive cancer clusters across patients using simplex plots.
- `11_Clinical_cox.ipynb`: Generates plots for the Cox proportional hazards model.

## Authors

- Lisa Fournier
- Garance Haefliger
- Albin Vernhes
- Lena Loye

## References

<a name="ref1"></a>
[1] Andersson, A., Larsson, L., Stenbeck, L., Salmén, F., Ehinger, A., Wu, S. Z., Al-Eryani, G., et al. (2021).
“Spatial Deconvolution of HER2-Positive Breast Cancer Delineates Tumor-Associated Cell Type Interactions.”
Nature Communications, 12(1), 6012.

<a name="ref2"></a>
[2] Thrane, K., Eriksson, H., Maaskola, J., Hansson, J., & Lundeberg, J. (2018).
“Spatially Resolved Transcriptomics Enables Dissection of Genetic Heterogeneity in Stage III Cutaneous Malignant Melanoma.”
Cancer Research, 78(20), 5970–5979.

<a name="ref3"></a>
[3] Ciga, O., Xu, T., & Martel, A. L. (2022).
“Self Supervised Contrastive Learning for Digital Histopathology.” 
Machine Learning with Applications, 7, 100198.

<a name="ref4"></a>
[4] Wang, X., Yang, S., Zhang, J., Wang, M., Zhang, J., Yang, W., Huang, J., & Han, X. (2022). 
“Transformer-Based Unsupervised Contrastive Learning for Histopathological Image Classification.”
Medical Image Analysis, 81, 102559.

<a name="ref5"></a>
[5] Chen, R. J., Ding, T., Lu, M. Y., Williamson, D. F. K., Jaume, G., Song, A. H., Chen, B., et al. (2024).
“Towards a General-Purpose Foundation Model for Computational Pathology.”
Nature Medicine, 30(3), 850–862.

<a name="ref6"></a>
[6] Xu, H., Usuyama, N., Bagga, J., Zhang, S., Rao, R., Naumann, T., Wong, C., et al. (2024).
“A Whole-Slide Foundation Model for Digital Pathology from Real-World Data.”
Nature. https://doi.org/10.1038/s41586-024-07441-w

<a name="ref7"></a>
[7] Zimmermann, E., et al. (2024).
“Virchow2: Scaling self-supervised mixed magnification models in pathology,”
arXiv preprint arXiv:2408.XXXXX.

<a name="ref8"></a>
[8] MahmoodLab. “UNI2-h.” GitHub repository. https://github.com/mahmoodlab/UNI

<a name="ref9"></a>
[9] Hörst, F., et al. (2024).
“CellViT: Vision Transformers for precise cell segmentation and classification.”
Medical Image Analysis, 94, 103143.

<a name="ref10"></a>
[10] Wang, X., et al. (2024).
“Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications.” 
Nature Communications, 15, 10232.
