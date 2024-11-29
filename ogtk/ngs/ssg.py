#!/usr/bin/env python3
import polars as pl
import argparse
import sys
from pathlib import Path
from typing import List, Tuple
from datetime import datetime
import rich

pprint = rich.print

def generate_all_combinations(df: pl.DataFrame) -> pl.DataFrame:
    """Generate all possible combinations of i5 and i7 indices using join."""
    # First filter for valid i5/i7 entries, then select ID
    i5_df = df.filter(pl.col('i5').is_not_null()).select('id')
    i7_df = df.filter(pl.col('i7').is_not_null()).select('id')
    
    # Cross join to get all combinations
    return (
        i5_df.join(i7_df, how='cross')
        .rename({'id': 'i5', 'id_right': 'i7'})
        .with_columns(pl.concat_str([pl.col('i5'), pl.col('i7')], separator='_').alias('sample'))
    )

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def process_indices(primer_db: Path, 
                   experiment_file: Path,
                   anchors: List[Tuple[str, str]],
                   project_name: str,
                   output_file: Path,
                   extract_regex: str= 'PO_(.*?)_',
                   i5_len: None|int= None,
                   i7_len: None|int= None,
                   i5_rc: bool=False,
                   generate_all: bool=False,
                   ) -> None|pl.DataFrame:
    """Process Illumina indices and generate sample sheet."""
    
    # Create anchors dataframe
    anchors_df = pl.DataFrame(
        anchors,
        schema=['upstream', 'downstream'],
        orient='row',
    )
    
    # Create regex patterns
    upstreams = "|".join(anchors_df.get_column("upstream"))
    downstreams = "|".join(anchors_df.get_column("downstream"))
    
    # Process primer database
    df = (
        pl.read_excel(primer_db)
        .with_columns(pl.col('id').str.extract(extract_regex, 1))
        .with_columns(
            pl.col('seq')
            .str.to_uppercase()
            .str.replace(f'^.*?({upstreams})', '')
            .alias('index')
        )
        .with_columns(
            pl.col('index')
            .str.to_uppercase()
            .str.replace(f'({downstreams}).*$', '')
        )
    )
    
    # Pivot the dataframe
    df = df.pivot(index='id', on='end', values='index')
    
    # If generate_all is True, create all combinations instead of reading experiment file
    pprint(":vampire:")

    if generate_all:
        exp = generate_all_combinations(df)
    else:
        exp = pl.read_excel(experiment_file)
    
    result = (
        exp
        .with_columns(pl.col('i5').cast(pl.Utf8))
        .with_columns(pl.col('i7').cast(pl.Utf8))
        .join(df.drop('i7'), left_on='i5', right_on='id')
        .join(df.drop('i5'), left_on='i7', right_on='id')
        .rename({'i5_right': 'index2', 'i7_right': 'index'})
        .drop('i5', 'i7')
    )
    
    result = (
            result
            .with_columns(pl.lit(project_name).str.replace_all(' ', '_').alias('Sample_Project'))
            .with_columns(pl.col('sample').str.replace_all(' ', '_'))
            )

    if i5_rc:
        result = result.with_columns(pl.col('index2').map_elements(reverse_complement, return_dtype=pl.Utf8))
    result = result.with_columns(pl.col('index').map_elements(reverse_complement, return_dtype=pl.Utf8))

    if i5_len is not None:
        result = result.with_columns(pl.col('index2').str.slice(0, i5_len))
    
    if i7_len is not None:
        result = result.with_columns(pl.col('index').str.slice(0, i7_len))

    # Select and rename columns for final output
    final = result.select([
        pl.col('sample').alias('Sample_ID'),
        pl.col('sample').alias('Sample_Name'),
        'index',
        'index2',
        'Sample_Project'
    ])
    
    include_full_header=False

    with open(output_file, 'w') as f:
        if include_full_header:
            f.write('[Header]\n')
            f.write('IEMFileVersion,5\n')
            f.write('Date,{}\n'.format(datetime.now().strftime('%Y-%m-%d')))

            f.write('Workflow,GenerateFASTQ\n')
            f.write('Application,FASTQ Only\n')
            f.write('Instrument Type,NovaSeq\n')
            f.write('Assay,TruSeq HT\n')
            f.write('Index Adapters,IDT-ILMN TruSeq DNA-RNA UD 96 Indexes\n')
            f.write('Chemistry,Amplicon\n\n')
            
            f.write('[Settings]\n')
            f.write('CreateFastqForIndexReads,1\n')
            f.write('BarcodeMismatchesIndex1,1\n')
            f.write('BarcodeMismatchesIndex2,1\n\n')
            
        f.write('[Data]\n')
        
        final.write_csv(f, include_header=True)
        return final

def validate_file(path: str) -> Path:
    """Validate that a file exists and return its Path object."""
    file_path = Path(path)
    if not file_path.exists():
        raise argparse.ArgumentTypeError(f"File {path} does not exist")
    return file_path

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate bcl2fastq sample sheet from primer and experiment data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--primer-db',
        type=validate_file,
        required=True,
        help='Path to primer database Excel file'
    )
    
    parser.add_argument(
        '--experiment',
        type=validate_file,
        help='Path to experiment Excel file'
    )
    
    parser.add_argument(
        '--project-name',
        type=str,
        required=True,
        help='Name of the sample project'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='SampleSheet.csv',
        help='Output path for sample sheet CSV'
    )
    
    parser.add_argument(
        '--i7_len',
        type=int,
        default=None,
        help='length of i7 index'
    )

    parser.add_argument(
        '--i5_len',
        type=int,
        default=None,
        help='length of i5 index'
    )
    
    parser.add_argument(
        '--i5_rc',
        action='store_true',
        help='Reverse complement i5'
    )

    parser.add_argument(
        '--all',
        action='store_true',
        help='Generate all combinations of indices. Ignores input experiment'
    )

    parser.add_argument(
        '--show',
        action='store_true',
        help='Show created sample sheet'
    )

    parser.add_argument(
        '--extr_regex',
        default='PO_(.*?)_',
        help='Regex to extract the first capture group: e.g. "PO_(.*?)_"'
    )

    # Add validation after parsing
    args = parser.parse_args()

    if not args.all and not args.experiment:
        parser.error("--experiment is required when not using --all")

    return args

def main():
    """Main function to run the sample sheet generator."""
    try:
        args = parse_args()
       
        anchors = [
            ('CGAGATCTACAC', 'ACACTCTTTC'),
            ('CGGCATACGAGAT', 'GTGACTGGAGTTC'),
        ]
        
        # Create Panel for arguments display
        from rich.panel import Panel
        from rich.text import Text
        from rich.align import Align

        text = Text()
        for key, value in vars(args).items():
            text.append(f"{key}: ", style="bold")
            text.append(f"{value}\n", style="cyan")
        
        rich.print(Panel(Align.center(text), title="[bold]Arguments[/bold]", border_style="blue"))



        final = process_indices(
            args.primer_db,
            args.experiment,
            anchors,
            args.project_name,
            Path(args.output),
            i5_len=args.i5_len,
            i7_len=args.i7_len,
            i5_rc=args.i5_rc,
            generate_all=args.all,
        )
        pprint(f"Successfully generated sample sheet: {args.output}")

        if args.show:
            print(final)

        
    except Exception as e:
        print(f"Error processing files: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
