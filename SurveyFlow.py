#!/usr/bin/env python3
"""
surveyFlow.py: A flexible pipeline for surveying genomic data using fastp, Jellyfish, GenomeScope, and Smudgeplot.
- Performs quality control, k-mer counting, genome feature estimation, and visualization.
- Uses .done files to track completed steps and skip them on subsequent runs.
Author: [HYG]
Date: March 31, 2025
"""
import subprocess
import os
import argparse
import logging
from datetime import datetime

def setup_logging(outdir, prefix):
    """设置日志，文件名基于 prefix 并输出到 outdir"""
    log_file = os.path.join(outdir, f"SurveyFlow.{prefix}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def run_command(cmd, step_name, done_file, logger):
    """执行命令并记录结果，成功后生成 .done 文件"""
    logger.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logger.info(f"Step {step_name} completed successfully.")
        # 成功后生成 .done 文件
        with open(done_file, 'w') as f:
            f.write(f"{step_name} completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode() if e.stderr else "No error message"
        logger.error(f"Error in {step_name}: {error_message}")
        raise

def fastp_quality_control(input_r1, input_r2, output_r1, output_r2, prefix, logger):
    """运行 fastp 进行质控，如果 .done 文件存在则跳过"""
    done_file = f"{prefix}_fastp.done"
    if os.path.exists(done_file):
        logger.info("Fastp step already completed (found .done file), skipping quality control step.")
        return
    cmd = (
        f"fastp -i {input_r1} -I {input_r2} -o {output_r1} -O {output_r2} -j {prefix}.fastp.json -h {prefix}.fastp.html"
    )
    run_command(cmd, "fastp quality control", done_file, logger)

def jellyfish_count(input_r1, input_r2, kmer_jf, threads, kmer_size, genome_size, logger):
    """运行 Jellyfish 进行 k-mer 计数，如果 .done 文件存在则跳过"""
    done_file = f"{kmer_jf}.done"
    if os.path.exists(done_file):
        logger.info("Jellyfish k-mer counting step already completed (found .done file), skipping counting step.")
        return
    cmd = (
        f"zcat {input_r1} {input_r2} | jellyfish count -C -m {kmer_size} -t {threads} -s {genome_size} -o {kmer_jf} /dev/stdin"
    )
    run_command(cmd, "Jellyfish k-mer counting", done_file, logger)

def jellyfish_histo(kmer_jf, kmer_histo, threads, logger):
    """生成 Jellyfish 直方图，如果 .done 文件存在则跳过"""
    done_file = f"{kmer_histo}.done"
    if os.path.exists(done_file):
        logger.info("Jellyfish histogram step already completed (found .done file), skipping histogram generation step.")
        return
    cmd = (
        f"jellyfish histo -t {threads} -o {kmer_histo} {kmer_jf}"
    )
    run_command(cmd, "Jellyfish histogram generation", done_file, logger)

def genomescope2(kmer_histo, prefix, threads, kmer_size, ploidy, logger):
    """运行 GenomeScope2 分析，如果 .done 文件存在则跳过"""
    output_dir = f"GenomeScope2_{prefix}"
    done_file = f"{output_dir}.done"
    if os.path.exists(done_file):
        logger.info("GenomeScope2 step already completed (found .done file), skipping analysis step.")
        return
    cmd = (
        f"genomescope2 -i {kmer_histo} -o {output_dir} -k {kmer_size} -m -1 -p {ploidy}"
    )
    run_command(cmd, "GenomeScope2 analysis", done_file, logger)

def smudgeplot_analysis(kmer_histo, kmer_jf, prefix, kmer_size, logger):
    """运行 Smudgeplot 分析，如果 .done 文件存在则跳过"""
    done_file = f"{prefix}_smudgeplot.done"
    if os.path.exists(done_file):
        logger.info("Smudgeplot step already completed (found .done file), skipping analysis step.")
        return
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
        run_command(hetkmers_cmd, "Smudgeplot hetkmers extraction", f"{prefix}_smudgeplot_hetkmers.done", logger)

        plot_cmd = (
            f"smudgeplot.py plot {prefix}_kmer_pairs_coverages.tsv -o {prefix}_smudgeplot -k {kmer_size}"
        )
        run_command(plot_cmd, "Smudgeplot visualization", done_file, logger)
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.decode() if e.stderr else "No error message"
        logger.error(f"Error in Smudgeplot analysis: {error_message}")
        raise

def check_software(logger):
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
    """创建输出目录，如果不存在则创建，并切换工作目录"""
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.chdir(outdir)  # 切换到输出目录作为工作目录

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

    # 设置输出目录并切换工作目录
    setup_output_directory(args.outdir)

    # 设置日志
    logger = setup_logging(args.outdir, args.prefix)

    # 检查软件依赖
    check_software(logger)

    # 文件名
    trimmed_r1 = f"{args.prefix}_trimmed_R1.fastq.gz"
    trimmed_r2 = f"{args.prefix}_trimmed_R2.fastq.gz"
    kmer_jf = f"{args.prefix}.jf"
    kmer_histo = f"{args.prefix}.histo"
    prefix = args.prefix

    # 运行流程
    try:
        fastp_quality_control(args.r1, args.r2, trimmed_r1, trimmed_r2, prefix, logger)
        jellyfish_count(trimmed_r1, trimmed_r2, kmer_jf, args.threads, args.kmer, args.size, logger)
        jellyfish_histo(kmer_jf, kmer_histo, args.threads, logger)
        genomescope2(kmer_histo, prefix, args.threads, args.kmer, args.ploidy, logger)
        smudgeplot_analysis(kmer_histo, kmer_jf, prefix, args.kmer, logger)
        logger.info("Genome survey pipeline completed successfully!")
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    main()
