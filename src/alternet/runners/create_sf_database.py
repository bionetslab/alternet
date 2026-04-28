from alternet.splicing_factors import *
import logging
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Integrate splicing factor data from SpliceAid-F, SplicingLore, KEGG, and ENCODE."
    )
    
    # Input Arguments
    parser.add_argument("--saf", required=True, help="Path to SpliceAid-F TSV file")
    parser.add_argument("--sl", required=True, help="Path to SplicingLore Excel file")
    parser.add_argument("--kegg", required=True, help="Path to KEGG JSON file")
    parser.add_argument("--encode", required=True, help="Path to ENCODE RBP Excel file")
    
    # Output Arguments
    parser.add_argument("-o", "--output", default="sf_database.csv", 
                        help="Path to save the resulting CSV (default: sf_database.csv)")

    args = parser.parse_args()

    try:
        # Build the DB
        db_df = build_sf_database(
            saf_path=args.saf,
            splicinglore_path=args.sl,
            kegg_path=args.kegg,
            encode_path=args.encode
        )

        # Save result
        db_df.to_csv(args.output, index=False)
        logging.info(f"Database successfully saved to: {args.output}")
        
    except Exception as e:
        logging.error(f"Failed to build database: {e}")

if __name__ == "__main__":
    main()