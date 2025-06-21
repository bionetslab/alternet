
<a id="readme-top"></a>





<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li><a href="#installation">Installation</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#important-files-and-folders">Important Files and Folders</a></li>


  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Gene regulatory networks (GRNs) play a role in understanding gene interactions and their regulatory effects. They help decode biological systems by identifying how genes interact and regulate cellular processes. However, conventional GRN inference methods operate at the gene-level, overlooking transcript-level variability introduced by Alternative-Splicing.
The goal of this pipeline is to infer gene regulatory networks that incorporate transcript-level information. This enables the discovery of regulatory mechanisms affected by alternative splicing, which are not captured in traditional gene-level analyses.

## Pipeline Workflow

<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="mi34qaba/SpliceAwareGRN/figures/MA_pipeline_page-0001.jpg" alt="Logo" width="80" height="80">
  </a>

The main pipeline script executes the following steps:

1. **Load Resources**  
   Reads Biomart annotations and a curated list of human transcription factors (TFs).

2. **Transcript Annotation Database**  
   Builds a transcript-level annotation database using APPRIS and DIGGER to support plausibility filtering.

3. **Data Preparation**  
   Processes expression data for both canonical (gene-level) and AS-aware (transcript-level) GRN inference.

   - `transcript_tfs`, `gene_tfs`, and `targets` must be Pandas DataFrames with:
     - **columns = samples**
     - **rows = Ensembl gene or transcript IDs**

4. **GRN Inference**  
   Infers GRNs using isoform-based and gene-based inputs. Multiple runs are aggregated for robustness.

5. **Isoform Categorization**  
   Classifies isoforms (e.g., dominant, balanced, non-dominant) based on their expression compared to other isoforms of the same gene.

6. **Aggregation and Filtering**  
   Aggregates GRNs from multiple runs and filters edges based on:
   - Feature importance
   - Edge frequency

7. **Edge Categorization**  
   Compares canonical and AS-aware networks to classify edges as:
   - Common
   - Gene-exclusive
   - Isoform-exclusive

8. **Plausibility Filtering**  
   Filters isoform-exclusive edges based on structural and functional annotations (APPRIS, DIGGER).


<p align="right">(<a href="#readme-top">back to top</a>)</p>



## Installation

#### Step 1: Clone the repo

   ```sh
   git clone https://github.com/HoffmannJuliane/SpliceAwareGRN.git
   ```

#### Step 2: Download Annotation files from APPRIS and DIGGER

* APPRIS
```bash
wget -P data/ https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/e110v48/appris_data.appris.txt
```

* DIGGER
```bash
wget -O data/digger_data.csv https://zenodo.org/records/3886642/files/domain_mapped_to_exons.csv

```

#### Step 3: Pixi Environment 

1. If not already, follow installation guide to install pixi:  https://pixi.sh/latest/installation/

2. Install environment:
  ```
  pixi install
  ```
3. Run in shell:
  ```
  pixi shell
  ```


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

### Step 1: Create Config File
  The config file saves all paths for the data needed for the transcript and annotation pipeline. A template of the config file is saved under /configs. With the script `create_config.py`you can create your own config based on the data you want to use.

  To create the config run:
  ```
  python create_config.py \
    --count_data path-to-gene-counts \
    --sample_attributes path-to-sample-attributes \
    --transcript_data path-to-transcript-counts \
    --parent_directory path-to-parent-directory \
    --results_dir /data/bionets/mi34qaba/SpliceAwareGRN/results/ \
    --nruns int
    --tissue tissue_name 
  ```

### Step 2:

Load Config File

```python

import yaml

with open(path-to-config-file, 'r') as f:
  config = yaml.safe_load(f)

```

### Step 3: Load Data 
The function for the inference and annotation pipeline requires three dataframes: `transcript_tfs`, `gene_tfs`, and `targets`.
`transcript_tfs`, `gene_tfs`, and `targets` must be Pandas DataFrames with:
  - **columns = samples**
  - **rows = Ensembl gene or transcript IDs**

### Step 4: Run Inference and Annotation Pipeline

```python
    from total_pipeline import *

    plausibility_filtered_df = inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets)
    
```

All networks are saved under the results_dir specified in the config file by default:

- aggregated AS-Aware GRN: SpliceAwareGRN/results/tissue-name/grn/tissue-name_as-aware.network_aggregated.tsv
- aggregated gene-level GRN: SpliceAwareGRN/results/tissue-name/grn/tissue-name_canonical.network_aggregated.tsv
- plausibility filtered isoform-unique GRN: SpliceAwareGRN/results/tissue-name/grn/tissue-name_as-aware.network_plausibility_filtered_iso_unique.tsv



<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Example Pipelines

### Inference with GTEX Data

If you are using GTex Data you can use the python script `total_pipeline.py` 

For example with Liver samples run in shell:
```
python total_pipeline.py -f ../configs/Liver.yaml
```



### Inference for MAGNet Data

The file `magnet.py`conducts the annotation and inference pipeline for MAGNet data and the file `downstream_analysis.py` shows an exemplary downstream analysis with a Gene Set Enrichment Analysis via gseapy.


## Important Files and Folders

- `/src/arboreto_added/`: Adjusted GRNBoost for transcript-level GRN Inference
- `create_config.py`: creates config yaml file needed as input for pipeline

- `total_pipeline.py`: shows total inference and annotation pipeline
- `inference.py`: isoform-level and gene-level GRN inference 
- `aggregate_grnboost.py`: functions to aggregate repeated GRNs to one consenus network
- `downstream_analysis.py`: Gene set enrichment analysis for plausibility filtered isoform-unique GRN based on gseapy
- `load_data.py`: functions to load GTEx expression data
- `utils_database.py`: helper functions for transcript annotation pipeline
- `utils_network.py`: helper functions for inference pipeline and GRN processing

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[product-screenshot]: images/screenshot.png

