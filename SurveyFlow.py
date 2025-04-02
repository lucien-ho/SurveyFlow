#!/usr/bin/env python3
"""
surveyFlow.py: A flexible pipeline for surveying genomic data using fastp, Jellyfish, GenomeScope, and Smudgeplot.
- Performs quality control, k-mer counting, genome feature estimation, and visualization.
Author: [HYG]
Date: March 31, 2025
"""
import subprocess
import os
import argparse
import logging
from datetime import datetime

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"surveyFlow_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_command(cmd, step_name):
    """执行命令并记录结果"""
    logger.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logger.info(f"Step {step_name} completed successfully.")
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode() if e.stderr else "No error message"
        logger.error(f"Error in {step_name}: {error_message}")
        raise

def fastp_quality_control(input_r1, input_r2, output_r1, output_r2, prefix):
    """运行 fastp 进行质控"""
    cmd = (
        f"fastp -i {input_r1} -I {input_r2} -o {output_r1} -O {output_r2} -j {prefix}.fastp.json -h {prefix}.fastp.html"
    )
    run_command(cmd, "fastp quality control")

def jellyfish_count(input_r1, input_r2, kmer_jf, threads, kmer_size, genome_size):
    """运行 Jellyfish 进行 k-mer 计数"""
    cmd = (
        f"zcat {input_r1} {input_r2} | jellyfish count -m {kmer_size} -t {threads} -s {genome_size} -o {kmer_jf} /dev/stdin"
    )
    run_command(cmd, "Jellyfish k-mer counting")

def jellyfish_histo(kmer_jf, kmer_histo, threads):
    """生成 Jellyfish 直方图"""
    cmd = (
        f"jellyfish histo -t {threads} -o {kmer_histo} {kmer_jf}"
    )
    run_command(cmd, "Jellyfish histogram generation")

def genomescope2(kmer_histo, prefix, threads, kmer_size, ploidy):
    """运行 GenomeScope2 分析"""
    cmd = (
        f"genomescope2 -i {kmer_histo} -o GenomeScope2_{prefix} -t {threads} -k {kmer_size} -m -1 -p {ploidy}"
    )
    run_command(cmd, "GenomeScope2 analysis")

def smudgeplot_analysis(kmer_histo, kmer_jf, prefix, kmer_size):
    """运行 Smudgeplot 分析，包括 cutoff、hetkmers 和 plot"""
    logger.info(f"Starting Smudgeplot analysis with {kmer_histo} and {kmer_jf}")
    try:
        l_cmd = f"smudgeplot.py cutoff {kmer_histo} L"
        l_result = subprocess.run(l_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        L = l_result.stdout.strip()
        logger.info(f"Lower cutoff (L): {L}")

        u_cmd = f"smudgeplot.py cutoff {kmer_histo} U"
        u_result = subprocess.run(u_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        U = u_result.stdout.strip()
        logger.info(f"Upper cutoff (U): {U}")

        hetkmers_cmd = (
            f"jellyfish dump -c -L {L} -U {U} {kmer_jf} | smudgeplot.py hetkmers -o {prefix}_kmer_pairs"
        )
        run_command(hetkmers_cmd, "Smudgeplot hetkmers extraction")

        plot_cmd = (
            f"smudgeplot.py plot {prefix}_kmer_pairs_coverages.tsv -o {prefix}_smudgeplot -k {kmer_size}"
        )
        run_command(plot_cmd, "Smudgeplot visualization")
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode() if e.stderr else "No error message"
        logger.error(f"Error in Smudgeplot analysis: {error_message}")
        raise

def check_software():
    """检查依赖软件是否可用"""
    required_software = {
        "fastp": "fastp --version",
        "jellyfish": "jellyfish --version",
        "genomescope2": "genomescope2 --version",
        "smudgeplot.py": "smudgeplot.py --version",
    }
    missing_tools = []

    for tool, command in required_software.items():
        try:
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                missing_tools.append(tool)
                logger.error(f"{tool} is not installed or not functioning correctly.")
        except FileNotFoundError:
            missing_tools.append(tool)
            logger.error(f"{tool} is not found in the system PATH.")

    if missing_tools:
        logger.error(f"Missing required tools: {', '.join(missing_tools)}. Please install them before running the pipeline.")
        exit(1)
    else:
        logger.info("All required tools are installed.")

def setup_output_directory(outdir):
    """创建输出目录，如果不存在则创建"""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        logger.info(f"Created output directory: {outdir}")
    return outdir

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="Genome analysis pipeline with fastp, Jellyfish, GenomeScope, and Smudgeplot.")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ file")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ file")
    parser.add_argument("--threads", type=int, default=16, help="Number of threads (default: 16)")
    parser.add_argument("--kmer", type=int, default=21, help="K-mer size (default: 21)")
    parser.add_argument("--prefix", default="output", help="Prefix for output files (default: output)")
    parser.add_argument("--ploidy", type=int, default=2, help="Ploidy (default: 2)")
    parser.add_argument("--size", required=True, help="Genome size")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
    args = parser.parse_args()

    check_software()

    # 设置输出目录并调整文件路径
    outdir = setup_output_directory(args.outdir)
    trimmed_r1 = os.path.join(outdir, f"{args.prefix}_trimmed_R1.fastq.gz")
    trimmed_r2 = os.path.join(outdir, f"{args.prefix}_trimmed_R2.fastq.gz")
    kmer_jf = os.path.join(outdir, f"{args.prefix}.jf")
    kmer_histo = os.path.join(outdir, f"{args.prefix}.histo")
    prefix_with_path = os.path.join(outdir, args.prefix)

    # 更新日志文件路径到输出目录
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.baseFilename = os.path.join(outdir, os.path.basename(handler.baseFilename))

    # 运行流程
    try:
        fastp_quality_control(args.r1, args.r2, trimmed_r1, trimmed_r2, prefix_with_path)
        jellyfish_count(trimmed_r1, trimmed_r2, kmer_jf, args.threads, args.kmer, args.size)
        jellyfish_histo(kmer_jf, kmer_histo, args.threads)
        genomescope2(kmer_histo, args.prefix, args.threads, args.kmer, args.ploidy)
        smudgeplot_analysis(kmer_histo, kmer_jf, prefix_with_path, args.kmer)
        logger.info("Genome survey pipeline completed successfully!")
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
