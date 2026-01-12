import sys
import os
import pytest
from unittest.mock import MagicMock, patch
import pandas as pd
import numpy as np
import tempfile

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from guide_prioritization import run_prioritization

@patch('guide_prioritization.predict_seq')
def test_prioritization_logic(mock_predict):
    # Setup real input file
    input_data = {
        'Spacer': ['A'*20, 'C'*20, 'G'*20],
        'targets': ['coding', 'noncoding', 'coding'],
        'GC': [0.5, 0.5, 0.9], 
        'tss_coord': [1000, 1000, 1000],
        'rs3_context': ['A'*30, 'C'*30, 'G'*30],
        'prom_start': [900, 900, 900],
        'prom_end': [1100, 1100, 1100],
        'offtarget_0_out_family': [0, 0, 5], 
        'gene': ['G1', 'G1', 'G1'],
        'polyT': [False, False, True] # Last one has polyT
    }
    df = pd.DataFrame(input_data)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_in:
        df.to_csv(tmp_in.name, sep='\t', index=False)
        tmp_in_path = tmp_in.name
        
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_out:
        tmp_out_path = tmp_out.name
    
    try:
        # Mock RS3 scores
        mock_predict.return_value = [10.0, 20.0, 30.0]
        
        run_prioritization(tmp_in_path, tmp_out_path)
        
        # Read output
        final_df = pd.read_csv(tmp_out_path, sep='\t')
        
        assert 'plant_act_score' in final_df.columns
        assert 'Tier' in final_df.columns
        assert 'is_offtarget_safe' in final_df.columns
        
        # Verify filtering didn't drop valid rows
        # Spacer A*20 (row 0), C*20 (row 1), G*20 (row 2).
        assert len(final_df) == 3
        
        # Logic Check
        r_polyC = final_df[final_df['Spacer'] == 'C'*20].iloc[0]
        assert "Tier 2" in r_polyC['Tier']
        assert r_polyC['is_offtarget_safe'] == True
        
        r_polyG = final_df[final_df['Spacer'] == 'G'*20].iloc[0]
        assert r_polyG['is_offtarget_safe'] == False
        # Should be Tier 5 (PolyT) because polyT=True overrides Tier 4 (unsafe) if it's top priority
        # Let's check my logic order: Tier 4 is assigned, then Tier 5. So it should be Tier 5.
        assert "Tier 5 (PolyT)" in r_polyG['Tier']
        
    finally:
        if os.path.exists(tmp_in_path): os.remove(tmp_in_path)
        if os.path.exists(tmp_out_path): os.remove(tmp_out_path)
