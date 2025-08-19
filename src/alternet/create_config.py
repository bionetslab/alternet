import argparse
import yaml
from pathlib import Path
import os.path as op

def main():
    parser = argparse.ArgumentParser(description="Create a YAML config file for SpliceAwareGRN.")

    parser.add_argument('--count_data', required=True, help='Path to GTEx gene TPM file')
    parser.add_argument('--sample_attributes', required=True, help='Path to GTEx sample attributes')
    parser.add_argument('--transcript_data', required=True, help='Path to GTEx transcript TPM file')
    parser.add_argument('--parent_directory', required=True, help='Path to parent directory')
    parser.add_argument('--results_dir', required=True, help='Path to results directory')
    parser.add_argument('--nruns', type=int, required=True, help='Number of inference repetitions')
    parser.add_argument('--tissue', required=True, help='Tissue name')
    

    args = parser.parse_args()

    config = {
        'biomart': '../data/biomart.txt',
        'count_data': args.count_data,
        'sample_attributes': args.sample_attributes,
        'transcript_data': args.transcript_data,
        'appris': '../data/appris_data.appris.txt',
        'digger': '../data/digger_data.csv',
        'parent_directory': args.parent_directory,
        'results_dir': args.results_dir,
        'tf_list': '../data/allTFs_hg38.txt',
        'nruns': args.nruns,
        'tissue': args.tissue
    }

    filename = op.join('../configs/', f'{args.tissue}.yaml')
    with open(filename, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"Config written to {'../configs/'}")

if __name__ == '__main__':
    main()

'''
example code 

python create_config.py \
  --count_data /data/bionets/mi34qaba/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct \
  --sample_attributes /data/bionets/mi34qaba/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt \
  --transcript_data /data/bionets/mi34qaba/data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct \
  --parent_directory /data/bionets/mi34qaba/ \
  --results_dir /data/bionets/mi34qaba/SpliceAwareGRN/results/ \
  --nruns 10 \
  --tissue Liver
'''