import argparse

import pandas as pd

from alternet.genome_annotation import read_annotation
from alternet.informative_exons import map_node_sources, node_exons, read_exon_nodes


def main():
    parser = argparse.ArgumentParser(
        description="Maps the Whippet node sources '<gene>:<node>' of a whippet_node_grn network to the "
                    "annotated exons that contain them."
    )
    parser.add_argument("--network", type=str, required=True, help="Network TSV from whippet_node_grn.")
    parser.add_argument("--gtf", type=str, required=True, help="GTF used for the Whippet index.")
    parser.add_argument("--exon_nodes", type=str, required=True,
                        help="Whippet exon-to-node mapping (graph.exons.tab(.gz) of the Whippet index).")
    parser.add_argument("--output", type=str, required=True, help="Output network TSV.")
    args = parser.parse_args()

    network = pd.read_csv(args.network, sep='\t')
    exons = node_exons(read_annotation(args.gtf, features={'exon'}), read_exon_nodes(args.exon_nodes))
    mapped = map_node_sources(network, exons)
    mapped.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(mapped)} rows for {len(network)} edges to {args.output} "
          f"({mapped['source_exon'].isna().sum()} rows without an annotated exon)")


if __name__ == "__main__":
    main()
