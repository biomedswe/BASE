from shortcuts import Shortcuts
from Misc2 import Misc


class cnv_calling:
    def __init__(self):
        self.shortcuts = Shortcuts()
        self.misc = Misc()

    def run_r_script(self, input_file_path, output_dir_path=None):
        if output_dir_path is None:
            output_dir_path = self.shortcuts.aligned_output_dir

        r_script_path = self.shortcuts.CNV_calling_path
        print(f"R script path: {r_script_path}")

        cmd = f"Rscript {r_script_path} {input_file_path} {output_dir_path}"
        self.misc.run_command(cmd)


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Process VCF files')
    parser.add_argument('-i', '--input', help='cnv reads path', required=True)
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    options = parser.parse_args()

    gvcf_parser = gvcf2cnv(options.input, options.output)
    gvcf_parser.run()


if __name__ == "__main__":
    main()
