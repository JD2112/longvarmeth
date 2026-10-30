import os
import sys
import csv
import json
import tempfile
import unittest
from pathlib import Path

# Ensure bin/ is on Python path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

from generate_methylation_dashboard import parse_methbat_tsv, calculate_stats


class TestSamplesheetIngestion(unittest.TestCase):
    """Validate samplesheet format, column requirements, and URI schemes."""

    REQUIRED_HEADERS = ["run_id", "input_path", "input_type", "kit", "model", "modified_bases_model"]
    VALID_INPUT_TYPES = ["pod5", "fast5", "bam", "fastq", "pacbio_bam", "pacbio_hifi", "hifi", "ubam", "aligned_bam"]

    def test_pacbio_samplesheet_structure(self):
        csv_path = REPO_ROOT / "tests" / "csv" / "samplesheet_pacbio_test.csv"
        self.assertTrue(csv_path.exists(), "PacBio test samplesheet must exist")

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            self.assertEqual(reader.fieldnames, self.REQUIRED_HEADERS)
            rows = list(reader)
            self.assertGreaterEqual(len(rows), 1)
            for row in rows:
                self.assertIn(row["input_type"], self.VALID_INPUT_TYPES)
                self.assertTrue(
                    row["input_path"].startswith("http://") or 
                    row["input_path"].startswith("https://") or 
                    Path(row["input_path"]).name.endswith(".bam")
                )

    def test_nanopore_samplesheet_structure(self):
        csv_path = REPO_ROOT / "tests" / "samplesheet_test.csv"
        self.assertTrue(csv_path.exists(), "Nanopore test samplesheet must exist")

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            self.assertEqual(reader.fieldnames, self.REQUIRED_HEADERS)
            rows = list(reader)
            self.assertGreaterEqual(len(rows), 1)
            for row in rows:
                self.assertIn(row["input_type"], self.VALID_INPUT_TYPES)


class TestMethylationDashboardGenerator(unittest.TestCase):
    """Validate parser, calculations, and HTML payload injection for the clinical dashboard."""

    def test_parse_methbat_tsv_with_mock_data(self):
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
            tf.write("chrom\tstart\tend\tcpg_label\texpected_category\tsummary_label\tmean_combined_methyl\tmean_hap1_methyl\tmean_hap2_methyl\tmean_meth_delta\tqc_warnings\n")
            tf.write("chr1\t1000\t2000\tH19-DMR\tImprinting\tPASS\t0.520\t0.980\t0.040\t0.940\tnone\n")
            tf.write("chr7\t5000\t6000\tPEG10\tImprinting\tPASS\t0.490\t0.950\t0.020\t0.930\tlow_coverage\n")
            tf.flush()
            temp_tsv = Path(tf.name)

        try:
            rows = parse_methbat_tsv(temp_tsv)
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["cpg_label"], "H19-DMR")
            self.assertEqual(rows[0]["qc_status"], 0)  # PASS
            self.assertEqual(rows[1]["qc_status"], 1)  # Warning
            self.assertEqual(float(rows[0]["mean_combined_methyl"]), 0.520)
            self.assertEqual(float(rows[0]["mean_hap1_methyl"]), 0.980)

            stats = calculate_stats(rows, "test_sample")
            self.assertEqual(stats["passCount"], 1)
            self.assertEqual(stats["cov"], "Phased-5mC")
            self.assertIn("96.5%", stats["h1Avg"])
        finally:
            if temp_tsv.exists():
                os.remove(temp_tsv)

    def test_parse_empty_tsv_fallback(self):
        with tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False) as tf:
            tf.write("chrom\tstart\tend\n")
            tf.flush()
            temp_tsv = Path(tf.name)

        try:
            rows = parse_methbat_tsv(temp_tsv)
            self.assertEqual(rows, [])
            stats = calculate_stats(rows, "empty_sample")
            self.assertEqual(stats["passCount"], 0)
            self.assertEqual(stats["cov"], "N/A")
        finally:
            if temp_tsv.exists():
                os.remove(temp_tsv)


class TestConfigAndModuleIntegrity(unittest.TestCase):
    """Validate nextflow configs, GIAB benchmarking profile, and subworkflow bindings."""

    def test_nextflow_profiles_exist(self):
        nf_config = REPO_ROOT / "nextflow.config"
        content = nf_config.read_text()
        for expected_profile in ["test_ont", "test_pacbio", "test_giab", "human", "eukaryote", "prokaryote", "virus", "fraka", "version_lock"]:
            self.assertIn(expected_profile, content, f"Profile '{expected_profile}' missing in nextflow.config")

    def test_assets_dashboard_template_present(self):
        template = REPO_ROOT / "assets" / "methylation_dashboard" / "dashboard-template.html"
        self.assertTrue(template.exists(), "dashboard-template.html must be present in assets")
        content = template.read_text(errors="ignore")
        self.assertIn("{{SAMPLE_DATA_JSON}}", content)
        self.assertIn("{{SAMPLE_STATS_JSON}}", content)

    def test_multiqc_config_contains_execution_status(self):
        mqc_config = REPO_ROOT / "assets" / "multiqc_config.yml"
        self.assertTrue(mqc_config.exists(), "multiqc_config.yml must exist")
        content = mqc_config.read_text()
        self.assertIn("pipeline_execution_status", content)


class TestTruthSetBenchmarkingIntegrity(unittest.TestCase):
    """Verify presence and validity of benchmark truth files and configuration."""

    def test_giab_truth_files_present(self):
        giab_dir = REPO_ROOT / "tests" / "data" / "giab"
        self.assertTrue(giab_dir.exists(), "GIAB data directory must exist")
        self.assertTrue((giab_dir / "hg002_truth_smallvar.vcf.gz").exists())
        self.assertTrue((giab_dir / "hg002_truth_smallvar.vcf.gz.tbi").exists())
        self.assertTrue((giab_dir / "hg002_benchmark_regions.bed").exists())
        self.assertTrue((giab_dir / "hg002_truth_sv.vcf.gz").exists())
        self.assertTrue((giab_dir / "hg002_truth_sv.vcf.gz.tbi").exists())
        self.assertTrue((giab_dir / "hg002_benchmark_sv_regions.bed").exists())

    def test_giab_config_properties(self):
        giab_cfg = REPO_ROOT / "conf" / "tests" / "test_giab.config"
        self.assertTrue(giab_cfg.exists(), "test_giab.config must exist")
        content = giab_cfg.read_text()
        self.assertIn("run_benchmarking", content)
        self.assertIn("benchmark_truth_vcf", content)
        self.assertIn("benchmark_sv_truth_vcf", content)


if __name__ == "__main__":
    unittest.main()
