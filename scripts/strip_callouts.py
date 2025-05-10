import re
import os
import shutil

def strip_callouts(input_path):
    """
    Strips callout blocks from a Quarto document and creates two files:
    - %filename%-report.qmd: Clean version without callouts
    - %filename%-teaching.qmd: Original file with callouts preserved
    """
    # Get the base path and filename without extension
    base_dir = os.path.dirname(input_path)
    filename = os.path.basename(input_path)
    filename_noext = os.path.splitext(filename)[0]
    
    # Define output paths
    teaching_path = os.path.join(base_dir, f"{filename_noext}-teaching.qmd")
    report_path = os.path.join(base_dir, f"{filename_noext}-report.qmd")
    
    # First, copy the original file to the teaching version
    shutil.copy2(input_path, teaching_path)
    
    # Now process the file to remove callouts
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

    print(f"✅ Original file copied to: {teaching_path}")
    print(f"✅ Stripped file saved as: {report_path}")

# Example usage
if __name__ == "__main__":
    input_qmd = "notebooks/veteran-lung-cancer-coxph.qmd"
    strip_callouts(input_qmd)
