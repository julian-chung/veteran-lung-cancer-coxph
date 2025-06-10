import re
import os
import shutil

def strip_callouts(input_path):
    """
    Strips callout blocks from a Quarto document.
    It creates two files:
    - %output_base%-teaching.qmd: Original file content with callouts preserved.
    - %output_base%-report.qmd: Clean version without callouts.
    %output_base% is derived from the input filename by removing a '-working' suffix if present.
    """
    # Get the base path and filename without extension
    base_dir = os.path.dirname(input_path)
    filename = os.path.basename(input_path)
    filename_noext = os.path.splitext(filename)[0]
    
    # Determine the base for output filenames
    # If filename_noext ends with "-working", remove it for the output base
    output_filename_base = filename_noext
    if filename_noext.endswith("-working"):
        output_filename_base = filename_noext[:-len("-working")]
    
    # Define output paths using the potentially modified base
    teaching_path = os.path.join(base_dir, f"{output_filename_base}-teaching.qmd")
    report_path = os.path.join(base_dir, f"{output_filename_base}-report.qmd")
    
    # First, copy the original file to the teaching version
    # The teaching file will have the new name but content from the original input_path
    shutil.copy2(input_path, teaching_path)
    
    # Now process the file to remove callouts
    # Read from the original input_path for stripping, not the teaching_path
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    output_lines = []
    inside_block = False

    # Updated pattern to detect start of callouts with .callout- syntax
    block_start_pattern = re.compile(r'^:::\s*\{\.callout-(note|tip|warning|important|caution)[^}]*\}')
    block_end_pattern = re.compile(r'^:::\s*$')

    for line in lines:
        if block_start_pattern.match(line):
            inside_block = True
            continue  # skip this line
        if inside_block and block_end_pattern.match(line):
            inside_block = False
            continue  # skip this line

        if not inside_block:
            output_lines.append(line)

    # Write the stripped content to the report file
    with open(report_path, 'w', encoding='utf-8') as f:
        f.writelines(output_lines)

    print(f"✅ Teaching version (original content) saved as: {teaching_path}")
    print(f"✅ Stripped report version saved as: {report_path}")

# Example usage
if __name__ == "__main__":
    # Adjust this to your target file name
    input_qmd = "notebooks/veteran-lung-cancer-coxph-working.qmd" 
    # Ensure this input file exists at the specified path relative to where you run the script
    if os.path.exists(input_qmd):
        strip_callouts(input_qmd)
    else:
        print(f"Error: Input file not found at {input_qmd}")
        print(f"Current working directory: {os.getcwd()}")