# Extended pre-training of histopathology foundation models uncovers co-existing breast cancer archetypes characterized by RNA splicing or TGF-β dysregulation

## Overview
Here, we develop a pipeline to systematically evaluate the biological concepts encoded within histopathology foundation models (hFMs) using molecular data. We also perform extended-pretraining of [UNI](https://github.com/mahmoodlab/UNI) to identify optimal conditions that enhance the model’s ability to encode richer, tumor tissue-specific biological concepts.

This code was developed on a 64-bit Debian GNU/Linux system (x86_64), running kernel version 6.1.0-26-amd64 (Debian 6.1.112-1), with dynamic preemption enabled (PREEMPT_DYNAMIC). The code to perform extended pretraining is available in the [dinov2](./dinov2/) folder, while the code to assess compute and evaluate the resulting embeddings is located in the [digitalhistopathology](./digitalhistopathology/) and [scripts](./scripts/) folders.

## 1. Installation

The environment can be installed via conda, with (`environment.yml`) or without (`environment_no_builds.yml`) platform-specific build constraints:

```bash
conda env create -f environment.yml
```

Please change the prefix before creating the environement at the end of the file.

**BE CAREFUL: If you want to run CTransPath you will a specific version of the timm library.**

Go to their [Github page](https://github.com/Xiyue-Wang/TransPath) under 1.CTransPath and install the modified timm library on top of the anaconda environment.
When you want to run other foundation models like Uni or ProvGigaPath, override the timm library with version 1.0.7:

```bash
pip install timm==1.0.7
```

If you want to be safer, create two anaconda environements starting from the yml file with one specific to CTransPath and one for the other models.

The typical installation time in around 25 minutes.

## 2. Data

### HER2-positive breast cancer dataset [[1](#ref1)]

The data can be downloaded from [Zenodo](https://zenodo.org/records/4751624).

- Rename the principal folder as `HER2_breast_cancer`
- Unzip all the folders: `count-matrices.zip`, `images.zip`, `meta.zip` and `spot-selections.zip`
- Place the principal folder in `data/`

All other new datasets can be added by creating a new class in `digitalhistopathology/datasets/real_datasets.py`


## 3. HuggingFace Hub

UNI and Prov-GigaPath are pretrained models whose access must be granted before using them. You must have a hugging face account, agree to the outlined terms of use, and request access.

Your hugging face access tokens (must remain private) should be set in `digitalhistopathology/access_token.py` or as read-only token as an environment variable: `export HF_TOKEN=<MY_READ_ONLY_TOKEN>`. If set, this value will overwrite the token stored on the machine (in either `HF_TOKEN_PATH` or `HF_HOME/token` if the former is not set). The first time loading the model, the weights will be downloaded to the hugging face hub cache in your home directory (`~/.cache/huggingface/hub`).

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

Here are the different steps of the pipeline that needs a model and a name to be initialized.

### 5.1. Compute patches HER2-positive breast cancer dataset:

First, we need to compute patches from the HER2 dataset. This can be done using the script `scripts/compute_patches_HER2.py`.
The resulting patches will be produced in the `results/compute_patches/her2_final_without_A/` folder. Patient A has been excluded from the analysis for suspicion of technical bias.

Example usage: 
```
python scripts/compute_patches_HER2.py
```

### 5.2. Compute image embeddings HER2-positive breast cancer dataset and k-NN selection of invasive cancer patches:

This step is made to compute the embeddings using a pre-trained model, and to apply k-NN on them to extract the patches classified as "invasive cancer patches". The shannon entropy computed on the explained variance distribution of SVD components is also computed, and molecular data are formatted for later use. This is done using the script `scripts/pipeline.py`. 

Here are the arguments you need to provide:
- `--model_name`: the pretrained model to be used. The supported values are `uni`, `ctranspath`, `provgigapath`, `simclr`, `virchow`. You can add other models in `digitalhistopathology/models.py`.
- `--retrained_model_path`: if you retrained the pre-trained model and want to extract the associated embeddings, specifiy the path to the weights file. 
- `--patches_folder`: folder to the patches dataset you want to compute the embeddings from. Default is `results/compute_patches/her2_final_without_A/`.
- `--pipeline_name`: how you want to name the pipeline
- `--results_folder`: folder to save the embeddings. 
All results will be saved at `{results_folder}/{pipeline_name}/`.
- `--dataset`: dataset name. Choose among "HER2" and "TNBC". Default is "HER2".


In more details here are the steps:
- Patches image embeddings will be computed and saved as `{results_folder}/{pipeline_name}/image_embedding.h5ad`:
- Select invasive cancer patches by computing KNN on `image_embeddding` by using label column as training data. At the end, we selected only invasive_cancer label to go more in depth into tumor heterogeneity.
   - Boxplot of the F1 and accuracy across patient to choose the optimal k of knn and barplots with patient and label fraction in each knn cluster are saved in the folder `{results_folder}/{pipeline_name}/select_invasive_cancer`.
   - The invasive images embedding is saved under `{results_folder}/{pipeline_name}/invasive_image_embedding.h5ad`.
- Shannon entropy is computed on the SVD decomposition and results are saved in the folder `{results_folder}/{pipeline_name}/shannon_entropy`.

Example usage to extract embeddings from SimCLR
```
python scripts/pipeline.py --model_name simclr --patches_folder ../results/compute_patches/her2_final_without_A/ --pipeline_name simclr --results_folder ../results/pipeline --embedding_name image_embedding 
```

Note: If you just want to extract embeddings without computing the shannon entropy or running the k-NN re-annotation, you can use the ready-to-use script `script/compute_embeddings.py` that takes the same arguments. 

This step requires a GPU. It was run on a rtx3090 GPU, 12 CPUs and 128GB of RAM. You can also compute without GPU as it is useful only for inference. It lasts around 20 minutes.

The config files for this step and all models are available under `config/pipeline/base_models`. You can refer to them to know the parameters you need to pass to the function.

### 5.3. Create the folder for invasive cancer patches only:

This step enables to create a new patches dataset, subset of the original dataset and filtered using a csv file. It is done with the script `script/create_subset_hdf5.py` that takes three arguments:

- `--original_hdf5`: Path to the original HDF5 file.
- `--csv_path`: Path to the CSV file containing the names of the patches to include.
- `--new_hdf5`: Path to save the new HDF5 file.


In the context of our study, we were interested in the invasive cancer patches, and after the kNN-relabeling of step 2, we created a filtered dataset to keep only the patches relabeled as "invasive cancer" from the pipeline applied to SimCLR embeddings.

```python
python scripts/create_subset_hdf5.py \
    --original_hdf5 ../results/compute_patches/her2_final_without_A/patches.hdf5 \
    --csv_path ../pipeline/simclr/select_invasive_cancer/selected_invasive_cancer.csv \
    --new_hdf5 ../results/compute_patches/filtered/patches.hdf5
``` 

This step is really fast and can be run on a single CPU with 16GB of memory. 

### 5.4. Extended-pretraining:

**5.4.1. Perform extended-pretraining**:

In this step, we performed extended pre-training of UNI, using the invasive cancer patches selected above. The final models are called T-UNI models. Different T-UNI models have been generated, varying the extended pretraining strategy (full or ExPLoRa), the loss (KDE or KoLeo), and finally, the number of prototypes. All used config files are available in `dinov2/configs`. The script to perform extended pre-training is `dinov2/dinov2/train/train.py` and can be used as follows:

```bash
PYTHONPATH=$(pwd) python dinov2/train/train.py path_to_your_config_file
```

For more information about the extended pretraining, please see the [dinov2 README](./dinov2/README.md). Note that you will need another environment to run it. The [model card](./dinov2/MODEL_CARD_T_UNI.md) is available for the T-UNI models. You can also download direclty the model weights on [Zenodo](https://zenodo.org/records/15053890) and place them in the folder `dinov2/extended_pretrained_models`. For the rest of the code to run, please put each model in a directory named with the model name and adding `HER2_` at the beginning. Then rename the weights file as `epoch_1.pth`.

Example: 
When you download the model `uni_full_kde_16384_prototypes.pth` on Zenodo, place it into `dinov2/extended_pretrained_models/HER2_uni_full_kde_16384_prototypes/epoch_1.pth`.

**5.4.2 Run the embedding extraction on the extended pre-trained models**:

Example usage to extract embeddings from a retrained version of UNI located in `dinov2/extended_pretrained_models/HER2_uni_full_kde_16384_prototypes` under `epoch_1.pth`:

```
python scripts/pipeline.py --model_name uni --retrained_model_path dinov2/extended_pretrained_models/HER2_uni_full_kde_16384_prototypes/epoch_1.pth --patches_folder results/compute_patches/her2_final_without_A/ --pipeline_name HER2_uni_full_kde_16384_prototypes --results_folder results/pipeline --embedding_name image_embedding
```

The config files to extract embeddings from all T-UNI models are available under `config/pipeline/T-UNI`.

### 5.5 Handcrafted features extraction

In this step, segmentation of the nuclei has been performed using CellViT [[9](#ref)]. Features related to the morphology and texture of the nuclei have been extracted using scMTOP57 and correspond to the “Nuclei-Morph”, and “Nuclei-Texture” features respectively. We have enhanced the “Nuclei-Morph” category by computing Zernike moments of each cell using the Mahotas python package58. We have adapted scMTOP to add color descriptors of the nuclei, referred to as “Nuclei-Color” features, namely mean, skewness, kurtosis and entropy of each RGB color channel, as well as intensity and transparency. Finally, statistics on the cell type composition have been computed based on cell labels outputted by CellViT, and are referred to as “Nuclei-Composition”. Features related to the extracellular matrix color (“ExtraCell-Color” features) and texture (“ExtraCell-Texture”) have been computed using scikit-image59. Similarly, texture and color features at the entire patch level have been extracted and are referred to as “WholePatch-Texture” and “WholePatch-Color”.

**5.5.1 Nuclei segmentation using CellViT**:

The [CellViT](./CellViT/) folder is a Git subtree cloned from the [original repository](https://github.com/TIO-IKIM/CellViT), with only minimal modifications made to adapt it to our environment. To segment the HER2 dataset using CellViT, you can use the script `digitalhistopathology/engineered_features/cell_segmentor.py` using the following arguments: 

To run the segmentation, you need to set up a specific environment for CellViT:

```bash
conda env create -f cellvit_env.yml
conda activate cellvit_env
```

Here are the key parameters for the `digitalhistopathology/engineered_features/cell_segmentor.py` script:
- `--segmentation_mode`: The segmentation mode to use. For CellViT, set this to `"cellvit"`.
- `--magnification`: The magnification level of the images (e.g., 20x).
- `--mpp`: Microns per pixel for the images.
- `--list_wsi_filenames`: List of paths to the wsi to process
- `--dataset_name`: Name of the dataset being processed.
- `--model_path`: Path to the pretrained CellViT model weights.
- `--results_saving_folder`: Directory where the segmentation results will be saved.

Usage:
```bash
python3 digitalhistopathology/engineered_features/cell_segmentor.py \
--segmentation_mode "cellvit" \
--magnification 20 \
--mpp 1 \
--list_wsi_filenames "../../data/HER2_breast_cancer/images/HE/B1.jpg ../../data/HER2_breast_cancer/images/HE/B2.jpg ../../data/HER2_breast_cancer/images/HE/B3.jpg ../../data/HER2_breast_cancer/images/HE/B4.jpg ../../data/HER2_breast_cancer/images/HE/B5.jpg ../../data/HER2_breast_cancer/images/HE/B6.jpg ../../data/HER2_breast_cancer/images/HE/C1.jpg ../../data/HER2_breast_cancer/images/HE/C2.jpg ../../data/HER2_breast_cancer/images/HE/C3.jpg ../../data/HER2_breast_cancer/images/HE/C4.jpg ../../data/HER2_breast_cancer/images/HE/C5.jpg ../../data/HER2_breast_cancer/images/HE/C6.jpg ../../data/HER2_breast_cancer/images/HE/D1.jpg ../../data/HER2_breast_cancer/images/HE/D2.jpg ../../data/HER2_breast_cancer/images/HE/D3.jpg ../../data/HER2_breast_cancer/images/HE/D4.jpg ../../data/HER2_breast_cancer/images/HE/D5.jpg ../../data/HER2_breast_cancer/images/HE/D6.jpg ../../data/HER2_breast_cancer/images/HE/E1.jpg ../../data/HER2_breast_cancer/images/HE/E2.jpg ../../data/HER2_breast_cancer/images/HE/E3.jpg ../../data/HER2_breast_cancer/images/HE/F1.jpg ../../data/HER2_breast_cancer/images/HE/F2.jpg ../../data/HER2_breast_cancer/images/HE/F3.jpg ../../data/HER2_breast_cancer/images/HE/G1.jpg ../../data/HER2_breast_cancer/images/HE/G2.jpg ../../data/HER2_breast_cancer/images/HE/G3.jpg ../../data/HER2_breast_cancer/images/HE/H1.jpg ../../data/HER2_breast_cancer/images/HE/H2.jpg ../../data/HER2_breast_cancer/images/HE/H3.jpg" \
--dataset_name "her2_final_without_A" \
--model_path "../../CellViT/models/pretrained/CellViT-SAM-H-x20.pth" \
--results_saving_folder "../../results/"
```
The results will be saved in `results/segmentation/CellViT/her2_final_without_A`.

This step took around 30 minutes on a h100 GPU with 12 cpus per task and 64GB of memory. 

**5.5.2 Handcrafted features computation**: 

This step computes handcrafted features such as morphology, texture, and color descriptors of nuclei and extracellular matrix. It also includes Zernike moments for nuclei and cell type composition statistics.

Here are the parameters for the `digitalhistopathology/engineered_features/engineered_features.py` script:
- `--method`: The feature extraction method to use (e.g., morphology, texture, color).
- `--path_to_cellvit_folder`: Path to the folder containing CellViT segmentation results.
- `--result_saving_folder`: Directory where the computed features will be saved.
- `--dataset_name`: Name of the dataset being processed.
- `--temporary_folder`: Temporary directory for intermediate computations.
- `--patches_info_filename`: File containing metadata about image patches.
- `--list_wsi_filenames`: List of whole slide image filenames to process.
- `--save_individual_wsi`: Flag to save results for each WSI individually.
- `--zernike_per_nuclei`: Flag to compute Zernike moments for each nucleus.
- `--num_cores`: Number of CPU cores to use for parallel processing.

Example usage:
```bash
python3 engineered_features.py \
--method scMTOP \
--path_to_cellvit_folder ../results/segmentation/CellViT/her2_final_without_A \
--result_saving_folder ../results/engineered_features/scMTOP/with_zernike \
--dataset_name her2_final_without_A \
--temporary_folder ../tmp \
--patches_info_filename ../results/computes_patches/her2_final_without_A/patches_info.pkl.gz \
--list_wsi_filenames "../../data/HER2_breast_cancer/images/HE/B1.jpg ../../data/HER2_breast_cancer/images/HE/B2.jpg ../../data/HER2_breast_cancer/images/HE/B3.jpg ../../data/HER2_breast_cancer/images/HE/B4.jpg ../../data/HER2_breast_cancer/images/HE/B5.jpg ../../data/HER2_breast_cancer/images/HE/B6.jpg ../../data/HER2_breast_cancer/images/HE/C1.jpg ../../data/HER2_breast_cancer/images/HE/C2.jpg ../../data/HER2_breast_cancer/images/HE/C3.jpg ../../data/HER2_breast_cancer/images/HE/C4.jpg ../../data/HER2_breast_cancer/images/HE/C5.jpg ../../data/HER2_breast_cancer/images/HE/C6.jpg ../../data/HER2_breast_cancer/images/HE/D1.jpg ../../data/HER2_breast_cancer/images/HE/D2.jpg ../../data/HER2_breast_cancer/images/HE/D3.jpg ../../data/HER2_breast_cancer/images/HE/D4.jpg ../../data/HER2_breast_cancer/images/HE/D5.jpg ../../data/HER2_breast_cancer/images/HE/D6.jpg ../../data/HER2_breast_cancer/images/HE/E1.jpg ../../data/HER2_breast_cancer/images/HE/E2.jpg ../../data/HER2_breast_cancer/images/HE/E3.jpg ../../data/HER2_breast_cancer/images/HE/F1.jpg ../../data/HER2_breast_cancer/images/HE/F2.jpg ../../data/HER2_breast_cancer/images/HE/F3.jpg ../../data/HER2_breast_cancer/images/HE/G1.jpg ../../data/HER2_breast_cancer/images/HE/G2.jpg ../../data/HER2_breast_cancer/images/HE/G3.jpg ../../data/HER2_breast_cancer/images/HE/H1.jpg ../../data/HER2_breast_cancer/images/HE/H2.jpg ../../data/HER2_breast_cancer/images/HE/H3.jpg" \--save_individual_wsi \
--zernike_per_nuclei \
--num_cores 12
```

The results will be saved in `results/engineered_features/scMTOP/with_zernike/her_final_without_A`. This step is particularly long to run, especially due to zernike moments computation. It can indeed take up to 7 hours per slide, with 16 CPUs and 24GB of memory. Therefore we recommand lauching the script in parallel for each WSI. In future versions, this script will include built-in parallelism to streamline the process, and the operations will be optimized in order to save computation time. 


### 5.6 Molecular data preparation

**5.6.1 Load molecular data**:
To load the molecular data properly, run the script:

```bash
python script/load_molecular_data.py --gene_embedding_saving_path ../results/molecular --patches_folder ../results/compute_patches/her2_final_without_A
```

This step is really fast (a few minutes) and can be run on a single cpu wih 16GB.

**5.6.2 Combat-correction**

In order to perform combat correction you need to run:

```bash
R scripts/combat_correction.R
```

This step is quite fast. It runs in 10minutes on a single cpu.

Note: The needeed R libraries are already installed in the conda environment installed at the beginning.

**5.6.2 Data formatting and UMAP computation of molecular embedding**:

In order to format the molecular embeddings and compute the UMAPs, run:

```bash
python scripts/molecular_formatting_and_UMAPs.py
```

This will create the `.h5ad`embeddings for the raw molecular data, the filtered molecular data, the filtered normalized molecular data and finally, the combat-corrected molecular data. 
This step is quite fast but requires a lot of memory. It runs in 10 minutes on a single cpu with 100GB of memory.

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

The notebooks analyze the results obtained in the full pipeline by doing some additional analyis or generating additional vizualisations. You can run them in their order of appearance. They will all load the needed files from a single config file. You can find it under `config/config_notebooks.json`.

**Description of the notebooks:**

- `1_Unsupervised_clustering_base.ipynb`: Analyze the results of the clustering and the ARI scores for the base models.
- `2_Regression_handcrafted_features.ipynb`: Analyze the results from linear regression, for both the base models and the retrained models. 
- `3_Unsupervised_clustering.ipynb`: Analyze the results of the clustering and the ARI scores for the retrained models. 
- `4_Invasive_cancer_clustering.ipynb`: Focuses on finding the best UMAP and k-means parameters for clustering invasive cancer patches to identify distinct archetypes. 
- `5_Invasive_clusters_image_wasserstein.ipynb`: Visualizes Wasserstein distances between invasive cancer clusters in the image embedding space. 
- `6_Invasive_clusters_molecular_wassertein.ipynb`: Visualizes Wasserstein distances between invasive cancer clusters in the molecular embedding space. 
- `7_Invasive_clusters_images.ipynb`: Visualizes images from invasive cancer clusters. 
- `8_Invasive_clusters_batch.ipynb`: Analyzes batch effects in invasive cancer clusters. 
- `9_Invasive_clusters_viz.ipynb`: Creates visualizations for invasive cancer clusters, on the original WSIs. 
- `10_Molecular_DGE_between_clusters.Rmd`: Performs differential gene expression analysis between identified molecular clusters. 
- `11_pathway_graphs.ipynb: Generates pathway enrichment graphs for molecular clusters. Note: prior to running this notebook, you need to go on [STRINGdb](https://string-db.org/) and put the markers genes generated from notebook 10. You need then to download the `.tsv` file related to the pathway you want to analyze. In this notebook, we show examples for 4 different pathways.   
- `12_gene_expression_viz.ipynb`: Visualizes gene expression patterns across clusters on the WSIs. 
- `13_Supp.ipynb`: Contains supplementary analysis showing why patient A has been removed from all the analysis. 


## Authors

- Lisa Fournier
- Garance Haefliger
- Albin Vernhes

## References

<a name="ref1"></a>
[1] Andersson, Alma, Ludvig Larsson, Linnea Stenbeck, Fredrik Salmén, Anna Ehinger, Sunny Z. Wu, Ghamdan Al-Eryani, et al. 2021. “Spatial Deconvolution of HER2-Positive Breast Cancer Delineates Tumor-Associated Cell Type Interactions.” Nature Communications 12 (1): 6012.

<a name="ref2"></a>
[2] Thrane, Kim, Hanna Eriksson, Jonas Maaskola, Johan Hansson, and Joakim Lundeberg. 2018. “Spatially Resolved Transcriptomics Enables Dissection of Genetic Heterogeneity in Stage III Cutaneous Malignant Melanoma.” Cancer Research 78 (20): 5970–79.

<a name="ref3"></a>
[3] Ciga, Ozan, Tony Xu, and Anne Louise Martel. 2022. “Self Supervised Contrastive Learning for Digital Histopathology.” Machine Learning with Applications 7 (March): 100198.

<a name="ref4"></a>
[4] Wang, Xiyue, Sen Yang, Jun Zhang, Minghui Wang, Jing Zhang, Wei Yang, Junzhou Huang, and Xiao Han. 2022. “Transformer-Based Unsupervised Contrastive Learning for Histopathological Image Classification.” Medical Image Analysis 81 (October): 102559.

<a name="ref5"></a>
[5] Chen, Richard J., Tong Ding, Ming Y. Lu, Drew F. K. Williamson, Guillaume Jaume, Andrew H. Song, Bowen Chen, et al. 2024. “Towards a General-Purpose Foundation Model for Computational Pathology.” Nature Medicine 30 (3): 850–62.

<a name="ref6"></a>
[6] Xu, Hanwen, Naoto Usuyama, Jaspreet Bagga, Sheng Zhang, Rajesh Rao, Tristan Naumann, Cliff Wong, et al. 2024. “A Whole-Slide Foundation Model for Digital Pathology from Real-World Data.” Nature, May. https://doi.org/10.1038/s41586-024-07441-w.

<a name="ref7"></a>
[7] Zimmermann, E., et al. “Virchow2: Scaling self-supervised mixed magnification models in pathology,” arXiv [cs.CV], 01-Aug-2024.

<a name="ref8"></a>
[8] MahmoodLab. “UNI2-h.” GitHub. https://github.com/mahmoodlab/UNI.

<a name="ref9"></a>
[9] Hörst, F., et al. “CellViT: Vision Transformers for precise cell segmentation and classification,” Med. Image Anal., vol. 94, p. 103143, May 2024.



