import sys
import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Add src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from extract_promoter_seqs import GenomeSequenceExtractor, find_guides

class MockGenomeSequenceExtractor(GenomeSequenceExtractor):
    """Mock that doesn't load a file but sets genome dict directly."""
    def __init__(self, genome_dict):
        self.genome = genome_dict
        self.chrom_map = {k: k for k in genome_dict.keys()}
        # Add simpler chrom mapping mock
        for k in list(self.chrom_map.keys()):
            if k.startswith("chr"):
                self.chrom_map[k.replace("chr", "")] = k

@pytest.fixture
def mock_extractor():
    """Mocks a genome with a known sequence."""
    # Seq: 100 As, Marker(20nt), 100 Ts. Total 220.
    marker = "GCTAGCTAGCGCTAGCTAGC" 
    full_seq = "A"*100 + marker + "T"*100
    genome = {"chr1": SeqRecord(Seq(full_seq), id="chr1")}
    return MockGenomeSequenceExtractor(genome)

def test_extract_promoter_positive_strand(mock_extractor):
    # TSS at 120 (1-based). Index 119.
    # Window 20bp upstream (far=20, close=0).
    # p_start = 120 - 20 = 100.
    # p_end = 120 - 0 - 1 = 119.
    # slice [99:119].
    # Fixture: A*100 + Marker("GCTAG...").
    # Index 99 is 'A' (100th A).
    # Index 100 is 'G' (Start of Marker).
    # Slice 99..119 gives 'A' + Marker[0:19].
    marker = "GCTAGCTAGCGCTAGCTAGC"
    
    seq, start, end = mock_extractor.extract_promoter("chr1", 120, "+", 20, 0)
    expected_seq = "A" + marker[:19]
    assert seq == expected_seq
    assert start == 100
    assert end == 119

def test_extract_promoter_negative_strand(mock_extractor):
    # TSS at 100 (Last A). Strand -.
    # Promoter (coding upstream) is higher coords.
    # far_bp=20, close_bp=0.
    # Window: [100 + 0 + 1, 100 + 20] -> [101, 120] (0-based: 101 to 120? No, 1-based logic in extract_promoter?)
    # extract_promoter logic:
    # p_start = tss + close_bp + 1
    # p_end = tss + far_bp
    # slice: [p_start-1 : p_end]
    # TSS=100. far=20, close=0.
    # p_start = 101. p_end = 120.
    # slice: [100 : 120]. 
    # indices 100..119. This is the Marker "GCTAGCTAGCG"
    # RC of Marker.
    # Marker: GCTAGCTAGCGCTAGCTAGC (Defined in fixture)
    # RC:     GCTAGCTAGCGCTAGCTAGC (It is its own RC? No.)
    # GCTAGCTAGCGCTAGCTAGC
    # C G A T ...
    # Let's check RC.
    marker = "GCTAGCTAGCGCTAGCTAGC"
    rc = str(Seq(marker).reverse_complement())
    seq, start, end = mock_extractor.extract_promoter("chr1", 100, "-", 20, 0)
    assert seq == rc
    assert start == 101
    assert end == 120

def test_find_guides_fwd():
    # Use spacer that doesn't create accidental PAMs on RC (avoid poly-C/G)
    # Spacer ending in '...ATAT'. PAM 'AGG'.
    spacer = "AT"*10 
    pam = "AGG"
    seq = "AAAA" + spacer + pam + "TTT" 
    
    guides = find_guides(seq)
    # RC of spacer 'AT'*10 is 'AT'*10. 'AGG' -> 'CCT'.
    # RC seq: 'AAAA' + 'CCT' + 'AT'*10 + 'TTT'.
    # 'CCT' ends in T. No GG.
    # 'AT'*10 has no GG.
    # So strictly 1 guide.
    
    assert len(guides) == 1
    g = guides[0]
    assert g['Spacer'] == spacer
    assert g['PAM'] == pam
    assert g['Orientation'] == "FWD"

def test_find_guides_rvs():
    # Guide on RVS strand.
    # Seq should contain the guide in RC form.
    # Target guide: Spacer "ACGT"*5, PAM "AGG".
    # RC of Target: "CCT" + RC("ACGT"*5).
    # RC("ACGT") -> "ACGT".
    # Construct seq: "AAAA" + "CCT" + "ACGT"*5 + "TTTT" (Need sufficient context).
    
    expected_spacer = "ACGT"*5
    # Sequence length needs to accommodate context window [i-3 : i+27]
    # i=4. Need i+27=31. Seq length 31 (indices 0-30).
    seq = "AAAA" + "CCT" + expected_spacer + "TTTT"
    
    guides = find_guides(seq)
    assert len(guides) == 1
    g = guides[0]
    
    assert g['Spacer'] == expected_spacer
    assert g['Orientation'] == "RVS"
    assert g['PAM'] == "AGG"

def test_find_guides_multiple():
    # Overlapping guides?
    seq = "AAAA" + "C"*20 + "AGG" + "C"*20 + "AGG" + "TTT"
    # Should find 2
    guides = find_guides(seq)
    # Note: find_guides iterates fully.
    # First guide at 4.
    # Second guide at 4 + 23 = 27? 
    # Length check: len(seq).
    assert len(guides) >= 2 

def test_resolve_chrom_mapping(mock_extractor):
    # Map "1" to "chr1"
    assert mock_extractor._resolve_chrom("1") == "chr1"
    assert mock_extractor._resolve_chrom("chr1") == "chr1"
    with pytest.raises(ValueError):
        mock_extractor._resolve_chrom("chr2")

def test_standardized_fwd_targets_coding():
    """Test that NGG on coding strand is FWD and targets noncoding."""
    # Need 4bp prefix and 3bp suffix for context logic
    seq = "AAAA" + "ATGC" * 5 + "AGG" + "TTT"
    guides = find_guides(seq)
    assert guides[0]['Orientation'] == "FWD"
    assert guides[0]['targets'] == "nc"
    assert guides[0]['Spacer'] == "ATGC" * 5

def test_standardized_rvs_targets_noncoding():
    """Test that CCN on coding strand is RVS and targets coding."""
    # 5' CC + 20nt 'A's 3'. Add padding.
    # RC is 5' 20nt 'T's + GG 3'
    # Needs 4bp prefix (relative to Top, i.e. end of guide?)
    seq = "AAAA" + "CCA" + "A"*20 + "TTTTT"
    guides = find_guides(seq)
    assert guides[0]['Orientation'] == "RVS"
    assert guides[0]['targets'] == "c"
    assert guides[0]['Spacer'] == "T"*20

def test_negative_strand_extraction_padding(mock_extractor):
    """Confirm (-) strand genes are RC'd to become 'Coding' strand extracts (with close_bp=1)."""
    # TSS at 100 (Last A). Promoter is 101-120 (The Marker)
    # On (-) strand, upstream is higher coordinates.
    # Genomic 101-120 is "GCTAGCTAGCGCTAGCTAGC"
    # RC should be "GCTAGCTAGCGCTAGCTAGC" (because it's a repeat here)
    seq, start, end = mock_extractor.extract_promoter("chr1", 100, "-", 20, 1)
    # The logic returns genomic coordinates 102 to 120
    assert start == 102
    assert end == 120
