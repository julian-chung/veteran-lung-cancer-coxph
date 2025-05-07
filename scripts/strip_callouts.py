import re

def strip_callouts(input_path, output_path):
    with open(input_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    output_lines = []
    inside_block = False

    # Pattern to detect start of callouts (tip, note, details, etc.)
    block_start_pattern = re.compile(r'^:::\s*\{\.?(tip|note|details)[^}]*\}')
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

    with open(output_path, 'w', encoding='utf-8') as f:
        f.writelines(output_lines)

    print(f"✅ Stripped file saved as: {output_path}")

# Example usage
if __name__ == "__main__":
    input_qmd = "notebooks/veteran-lung-cancer-coxph.qmd"
    output_qmd = input_qmd.replace(".qmd", "_report.qmd")
    strip_callouts(input_qmd, output_qmd)
