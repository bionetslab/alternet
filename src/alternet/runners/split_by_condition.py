import argparse
import os
import os.path as op
import re
import shutil
from glob import glob

import pandas as pd


def condition_dir_name(condition):
    """Directory name for a condition: characters other than letters, digits, '.', '_' and '-' become '_'."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', str(condition))


def read_mapping(mapping_file, sample_column=None, condition_column=None):
    """
    Reads the sample -> condition mapping.

    Parameters:
        mapping_file (str): TSV with a header, one row per sample.
        sample_column (str): Column with the run IDs (e.g. SRR...). First column if None.
        condition_column (str): Column with the conditions. Second column if None.

    Returns:
        pd.Series: condition per sample (index), without rows that miss either value.
    """
    mapping = pd.read_csv(mapping_file, sep='\t', dtype=str)
    sample_column = sample_column or mapping.columns[0]
    condition_column = condition_column or mapping.columns[1]
    mapping = mapping.dropna(subset=[sample_column, condition_column])
    duplicated = mapping.loc[mapping[sample_column].duplicated(), sample_column]
    if len(duplicated):
        raise ValueError(f"Samples with several rows in {mapping_file}: {sorted(set(duplicated))}")
    return mapping.set_index(sample_column)[condition_column]


def split_by_condition(data_dir, mapping, output_dir):
    """
    Splits a gene expression folder by condition into one folder per condition with the same layout.

    The input folder has one directory per sample (<sample>/counts/...) and merged tables with one
    column per sample in counts/ (*.csv, *.tsv). For each condition, output_dir/<condition>/ gets
    - <sample>/: the sample directory, moved out of data_dir, for each of the condition's samples,
    - counts/: the merged tables with only the condition's sample columns (same file names).
    Non-numeric columns of the merged tables (e.g. 'Name') are kept in every condition's table.
    The merged tables in data_dir/counts/, the sample directories without a condition and other files
    and folders at the top level of data_dir (e.g. logs/, qc/) stay in data_dir. Nothing is moved if a
    sample directory already exists in output_dir.

    Parameters:
        data_dir (str): Gene expression folder (e.g. GSE130036/).
        mapping (pd.Series): Condition per sample (see read_mapping).
        output_dir (str): Output folder, one subfolder per condition.

    Returns:
        pd.DataFrame: One row per mapped sample with columns 'sample', 'condition', 'sample_dir'
            (whether it has a sample directory) and, per merged table, whether it is a column of it.
    """
    sample_dirs = {op.basename(d) for d in glob(op.join(data_dir, '*')) if op.isdir(op.join(d, 'counts'))}
    unmapped = sorted(sample_dirs - set(mapping.index))
    if unmapped:
        print(f"{len(unmapped)} sample directories without a condition, not moved: {unmapped}")
    other = sorted(set(os.listdir(data_dir)) - sample_dirs - {'counts'})
    if other:
        print(f"Not moved from the top level of {data_dir}: {other}")

    report = pd.DataFrame({'sample': mapping.index, 'condition': mapping.values})
    report['sample_dir'] = report['sample'].isin(sample_dirs)
    moves = {op.join(data_dir, sample): op.join(output_dir, condition_dir_name(condition), sample)
             for sample, condition in mapping.items() if sample in sample_dirs}
    existing = [dest for dest in moves.values() if op.exists(dest)]
    if existing:
        raise FileExistsError(f"{len(existing)} sample directories already exist in {output_dir}, e.g. {existing[0]}")

    tables = sorted(glob(op.join(data_dir, 'counts', '*.csv')) + glob(op.join(data_dir, 'counts', '*.tsv')))
    for table_file in tables:
        sep = '\t' if table_file.endswith('.tsv') else ','
        table = pd.read_csv(table_file, sep=sep)
        annotation = list(table.columns.difference(table.select_dtypes('number').columns, sort=False))
        report[op.basename(table_file)] = report['sample'].isin(table.columns)
        for condition, group in mapping.groupby(mapping):
            out_dir = op.join(output_dir, condition_dir_name(condition), 'counts')
            os.makedirs(out_dir, exist_ok=True)
            columns = annotation + [s for s in group.index if s in table.columns]
            table[columns].to_csv(op.join(out_dir, op.basename(table_file)), sep=sep, index=False)

    for source, dest in moves.items():
        os.makedirs(op.dirname(dest), exist_ok=True)
        shutil.move(source, dest)

    return report


def main():
    parser = argparse.ArgumentParser(
        description="Splits a gene expression folder (<sample>/ directories and merged tables in counts/) by "
                    "condition. Per condition: its <sample>/ directories (moved) and the merged tables with "
                    "only its sample columns."
    )
    parser.add_argument("--data_dir", type=str, required=True, help="Gene expression folder, e.g. GSE130036/.")
    parser.add_argument("--mapping", type=str, required=True,
                        help="TSV with a header and one row per sample: sample ID (e.g. SRR...) and condition.")
    parser.add_argument("--sample_column", type=str, required=False,
                        help="Column of --mapping with the sample IDs. Default: first column.")
    parser.add_argument("--condition_column", type=str, required=False,
                        help="Column of --mapping with the conditions. Default: second column.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output folder, one subfolder per condition.")
    args = parser.parse_args()

    mapping = read_mapping(args.mapping, args.sample_column, args.condition_column)
    report = split_by_condition(args.data_dir, mapping, args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)
    report.to_csv(op.join(args.output_dir, 'samples.tsv'), sep='\t', index=False)

    print(report.groupby('condition').size().rename('n_samples').to_string())
    missing = report.loc[~report.drop(columns=['sample', 'condition']).all(axis=1)]
    if len(missing):
        print(f"{len(missing)} mapped samples are missing from at least one table or have no sample directory "
              f"(see {op.join(args.output_dir, 'samples.tsv')}):")
        print(missing.to_string(index=False))


if __name__ == "__main__":
    main()
