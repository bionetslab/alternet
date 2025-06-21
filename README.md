
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
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#important-files-and-folders">Important Files and Folders</a></li>


  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

AS-Aware Gene Regulatory Network inference


<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Installation

#### Step 1: Clone 

Clone the repo
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

### Step 3:

Run inference and annotation pipeline by calling:

```python

    transcript_tfs, gene_tfs, targets = load_data 
    
    inference_and_annotation_pipeline(config, transcript_tfs, gene_tfs, targets)
    
```




<p align="right">(<a href="#readme-top">back to top</a>)</p>


## Important Files and Folders

- `/src/arboreto_added/`: Adjusted GRNBoost for transcript-level GRN Inference
- `create_config.py`: creates config yaml file needed as input for pipeline

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[product-screenshot]: images/screenshot.png

