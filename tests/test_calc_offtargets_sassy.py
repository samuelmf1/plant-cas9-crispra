import sys
import os
import pytest
from unittest.mock import MagicMock, patch
import pandas as pd
from Bio.Seq import Seq

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from calc_offtargets_sassy import get_trimmed_seq, parse_chrom_from_text_id, check_overlap, run_sassy_and_count, load_family_regions

def test_get_trimmed_seq():
    # Helper to build 30bp seq: 4pre + 23spacer + 3suf
    prefix = "AAAA"
    suffix = "TTT"
    core = "C"*20 + "AGG" # 23bp
    seq = prefix + core + suffix
    
    # 1. (+, FWD) -> core is Top. No RC.
    assert get_trimmed_seq(seq, "FWD", "+") == core
    
    # 2. (+, RVS) -> core is Bottom. Need RC.
    assert get_trimmed_seq(seq, "RVS", "+") == str(Seq(core).reverse_complement())
    
    # 3. (-, FWD) -> core is Bottom. Need RC.
    assert get_trimmed_seq(seq, "FWD", "-") == str(Seq(core).reverse_complement())
    
    # 4. (-, RVS) -> core is Top. No RC.
    assert get_trimmed_seq(seq, "RVS", "-") == core
    
    # Short seq (should return as is)
    short = "ACGT"
    assert get_trimmed_seq(short, "FWD", "+") == short

def test_parse_chrom_from_text_id():
    # Code expects Ensembl/TAIR format: chromosome:TAIR10:1:1:30427671:1
    # regex: chromosome:[^:]+:([^:]+):
    valid_id = "dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF"
    assert parse_chrom_from_text_id(valid_id) == "1"
    
    # Fallback
    assert parse_chrom_from_text_id("chr1") == "chr1"

def test_check_overlap():
    # Region: (chrom, start, end)
    # check_overlap(chrom, s, e, regions)
    regions = [("chr1", 100, 200), ("chr2", 500, 600)]
    
    # Same chrom, overlap
    assert check_overlap("chr1", 150, 160, regions) == True
    assert check_overlap("chr1", 90, 110, regions) == True
    assert check_overlap("chr1", 190, 210, regions) == True
    assert check_overlap("chr1", 100, 200, regions) == True
    
    # Same chrom, no overlap
    assert check_overlap("chr1", 10, 99, regions) == False
    assert check_overlap("chr1", 201, 300, regions) == False
    
    # Diff chrom
    assert check_overlap("chr2", 150, 160, regions) == False
    assert check_overlap("chr2", 550, 560, regions) == True

def test_load_family_regions():
    df = pd.DataFrame({
        'family_id': ['A', 'A', 'B'],
        'chrom': ['chr1', 'chr1', 'chr2'],
        'start': [100, 500, 100],
        'end': [200, 600, 200]
    })
    fam_map = load_family_regions(df)
    assert len(fam_map['A']) == 2
    assert len(fam_map['B']) == 1
    assert ('chr1', 100, 200) in fam_map['A']

@patch('subprocess.Popen')
def test_run_sassy_and_count(mock_popen):
    # Mock sassy output:
    # pat_id text_id cost ... start end
    # map pat_id to seq.
    
    # Correct columns? Sassy output format varies.
    # Code dynamically maps columns by header.
    header = "pat_id\ttext_id\tcost\tmatches\tstart\tend\tstrand\n"
    # pat_id=0 (matches input id 0). cost=0.
    line1 = "0\tchr1:1000-2000(+)|gene:G1\t0\t23\t1050\t1073\t+\n"
    # pat_id=1. cost=1.
    line2 = "1\tchr1:1000-2000(+)|gene:G1\t1\t22\t1080\t1103\t+\n"
    
    process_mock = MagicMock()
    process_mock.stdout.__iter__.return_value = [header, line1, line2]
    process_mock.returncode = 0
    process_mock.stderr.read.return_value = ""
    mock_popen.return_value = process_mock
    
    pat_to_indices = {"SEQ_A": [0], "SEQ_B": [1]}
    id_to_pat = {0: "SEQ_A", 1: "SEQ_B"}
    
    df_fams = pd.Series([None, None], index=[0, 1]) # No family logic
    
    counts_tot, c_in, c_out = run_sassy_and_count(
        ["sassy"], pat_to_indices, id_to_pat, df_fams, {}, False
    )
    
    # SEQ_A (idx 0) has 1 hit at cost 0.
    assert counts_tot[0][0] == 1
    # SEQ_B (idx 1) has 1 hit at cost 1.
    assert counts_tot[1][1] == 1

@patch('subprocess.Popen')
def test_offtarget_0_self_hit_check(mock_popen):
    # Verify that a pattern matching the genome with 0 mismatch_intolerance yields offtarget_0 >= 1
    # User Requirement: every guide should have at least 1 offtarget_0 (self hit).
    
    header = "pat_id\ttext_id\tcost\tmatches\tstart\tend\tstrand\\n"
    # Pattern 1 matches itself at cost 0 (Self-hit)
    line1 = "1\tchr1:100\t0\t23\t100\t123\t+\\n"
    # Pattern 1 also matches elsewhere at cost 0 (Duplicate?)
    line2 = "1\tchr1:500\t0\t23\t500\t523\t+\\n"
    
    process_mock = MagicMock()
    process_mock.stdout.__iter__.return_value = [header, line1, line2]
    process_mock.returncode = 0
    process_mock.stderr.read.return_value = ""
    mock_popen.return_value = process_mock
    
    pat_to_indices = {"ACGT": [10]} # Row 10 has this seq
    id_to_pat = {1: "ACGT"} # sassy ID 1 maps to this seq
    
    counts_tot, _, _ = run_sassy_and_count(["sassy"], pat_to_indices, id_to_pat, {}, {}, False)
    
    # Row 10 should have 2 hits at cost 0
    assert counts_tot[10][0] == 2
    # Requirement: At least 1
    assert counts_tot[10][0] >= 1
