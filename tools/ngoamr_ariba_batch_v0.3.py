import os
import re
import argparse
import subprocess

## ngoamr_ariba_batch_v0.3.0.py (by Leonor Sánchez Busó, last updated 06/05/2026)
## Script to run ARIBA on .fastq.gz files across one or multiple directories
## Customized to work with a custom ARIBA database for N. gonorrhoeae
## The extract_indels() function creates a modified report.tsv file reporting the 
## penA insD345 insertion and mtrR_promoter -53Adel deletions so they can be called by 'ariba summary'
## Execute as: python ~/github/NgoAMR/ngoamr_ariba_batch_v0.3.py -d /path1/to/fastq_files,/path2/to/fastq_files -o <OUTPUT_DIR> --db_path ARIBA_NgoDB_040526

## Functions ##

def execute_ariba_run(threads, db_path, forward_file, reverse_file, output_folder):
    # Run the ARIBA command with the specified number of threads
    cmd = ['ariba', 'run', '--verbose', '--threads', str(threads), db_path, forward_file, reverse_file, output_folder]
    subprocess.run(cmd)

def execute_ariba_summary(output_dir, short_names):
    print(f"# Running ARIBA summary...")
    # Generate filenames.txt for the ARIBA summary including absolute paths for compatibility with sensityping
    with open('filenames.txt', 'w') as f:
        report_files = [os.path.join(output_dir, d, 'report_complete.tsv') for d in os.listdir(output_dir) if d.endswith('_ARIBA')]
        if short_names:
            report_names = [d.replace('_ARIBA', '') for d in os.listdir(output_dir) if d.endswith('_ARIBA')]
            for count, report_file in enumerate(report_files):
                f.write(report_file+'\t'+report_names[count]+'\n')
        else:
            for count, report_file in enumerate(report_files):
                f.write(report_file+'\n')
    # Execute ARIBA summary
    cmd = [
        'ariba', 'summary', 'ariba_summary', '-f', 'filenames.txt', '--cluster_cols',
        'assembled,ref_seq,pct_id', '--col_filter', 'n',
        '--row_filter', 'n', '--no_tree', '--v_groups', '--known_variants'
    ]
    subprocess.run(cmd)

def extract_indels(output_folder):
    # Modify the report.tsv file to include the penA.insD345 and mtrR promoter deletion if detected
    report_file = os.path.join(output_folder, 'report.tsv')
    new_report_file = os.path.join(output_folder, 'report_complete.tsv')
    with open(report_file, 'r') as infile, open(new_report_file, 'w') as outfile:
        for line in infile:
            if 'D147_T148insT' in line:
                pattern = r"0\t\.\tp\t\.\t0\tD147_T148insT"
                replacement = "1\tSNP\tp\tD147_T148insT\t1\tD147_T148insT"
                modified_line = re.sub(pattern, replacement, line)
            elif 'R146_D147insR' in line:
                pattern = r"0\t\.\tp\t\.\t0\tR146_D147insR"
                replacement = "1\tSNP\tp\tR146_D147insR\t1\tR146_D147insR"
                modified_line = re.sub(pattern, replacement, line)
            elif 'A197.' in line: #0    .   n   .   0   A197.
                pattern = r"0\t\.\tn\t\.\t0\tA197\."
                replacement = "1\tSNP\tn\tA197.\t1\tA197."
                modified_line = re.sub(pattern, replacement, line)
            else:
                modified_line = line
            outfile.write(modified_line)

def run_wrapper(input_dirs, output_dir, db_path, threads):
    # Split the input directories by comma
    dirs = input_dirs.split(',')

    # Define possible forward and reverse suffix patterns
    forward_suffixes = ['_1.fastq.gz', '_R1.fastq.gz', '_R1_001.fastq.gz', '_forward_paired.fq.gz']
    reverse_suffixes = ['_2.fastq.gz', '_R2.fastq.gz', '_R2_001.fastq.gz', '_reverse_paired.fq.gz']

    # Loop through each directory
    for d in dirs:
        # Get all the forward read files in the directory based on possible suffixes
        forward_files = [f for f in os.listdir(d) for suffix in forward_suffixes if f.endswith(suffix)]

        for f in forward_files:
            # Extract the sample name based on the detected suffix
            for suffix in forward_suffixes:
                if f.endswith(suffix):
                    sample_name = f.rsplit(suffix, 1)[0]
                    reverse_suffix = reverse_suffixes[forward_suffixes.index(suffix)]
                    print(f"# Running {sample_name}...")
                    break

            # Construct the reverse read file name
            reverse_file = os.path.join(d, sample_name + reverse_suffix)

            # Construct each output folder inside the main output directory
            output_folder = os.path.join(output_dir, sample_name + '_ARIBA')

            # Check if the output folder and report.tsv already exists
            # If it exists and the report.tsv file does not exist, delete the folder so ARIBA can be re-run without errors 
            if os.path.exists(output_folder):
                if os.path.exists(output_folder+'/report_complete.tsv'):
                    print(f"#    {sample_name} completed.")
                else:
                    print(f"#    {sample_name} already exists but report_complete.tsv is not found")
                    if os.path.exists(output_folder+'/report.tsv'):
                        if args.modify_report:
                            print(f"#    report.tsv is present, extracting key indels into report_complete.tsv...")
                            # Tune the report.tsv file so the penA insD345 and mtrR -53Adel determinants are reported later by ARIBA summary
                            extract_indels(output_folder)
                    else:
                        print(f"#    Deleting existing folder and re-running ARIBA in this sample...")
                        del_cmd = [ 'rm', '-r', output_folder]
                        subprocess.run(del_cmd)
                        execute_ariba_run(threads, db_path, os.path.join(d, f), reverse_file, output_folder)
                        # Tune the report.tsv file so the penA insD345 and mtrR -53Adel determinants are reported later by ARIBA summary
                        if not os.path.exists(output_folder+'/report_complete.tsv'):
                            extract_indels(output_folder)
            else:
                execute_ariba_run(threads, db_path, os.path.join(d, f), reverse_file, output_folder)
                # Tune the report.tsv file so the penA insD345 and mtrR -53Adel determinants are reported later by ARIBA summary
                if not os.path.exists(output_folder+'/report_complete.tsv'):
                    extract_indels(output_folder)
            
            # Run ariba summary
            execute_ariba_summary(output_dir, args.short_names)

    
## Main ##

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NgoAMR ARIBA_batch: run ARIBA on multiple .fastq.gz files and directories")
    parser.add_argument('-d', '--input_dirs', required=True, help="Comma-separated list of input directories containing .fastq.gz files.")
    parser.add_argument('-p', '--db_path', required=True, help="Path to the ARIBA database.")
    parser.add_argument('-m', '--modify_report', required=False, help="Only create a report_complete.tsv file with key indels extracted from existing report.tsv.", action='store_true')
    parser.add_argument('-s', '--summary_only', required=False, help="Run ariba summary only.", action='store_true')
    parser.add_argument('-n', '--short_names', required=False, help="Use folder names and not full paths in the summary report. Note: full paths are needed for sensityping (separate tool).", action='store_true')
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for ARIBA results.")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads for ARIBA run.")

    args = parser.parse_args()

    # Construct the output folder
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        create_folder = ['mkdir', output_dir]
        subprocess.run(create_folder)
    output_dir = os.path.abspath(output_dir)
    print(output_dir)

    if args.summary_only:
        execute_ariba_summary(output_dir, args.short_names)
    else:
        run_wrapper(args.input_dirs, output_dir, args.db_path, args.threads)

