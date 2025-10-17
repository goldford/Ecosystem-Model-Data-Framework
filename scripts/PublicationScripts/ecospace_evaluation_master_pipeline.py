"""
Script:   Master evaluation pipeline for ecospacec outputs analysis
Created by: G Oldford, 2025
Purpose:
   - orchestrates running of evaluation by calling functions etc in other scripts
Input:
  - config file
  - calls several eval scripts
Output:
  - log files
  - results of others scripts
"""


import ecospace_eval_1a_ASC_to_NC as data_prep
import ecospace_eval_2a_mcewan_phyto_seas as phyto_seas
import ecospace_eval_2b_nemcek_phyto_match as phyto_match_1
import ecospace_eval_2c_nemcek_vs_ecospace as phyto_eval_1
import ecospace_eval_3_QU39_match as phyto_match_2
import ecospace_eval_3b_QU39_vs_ecospace as phyto_eval_2
import ecospace_eval_4a_bloomt_CSoG as bt_eval_1
import ecospace_eval_6_nutrients as nu_eval
import ecospace_eval_5_zoop_eval as zp_eval



def run_ecospace_pipeline(run_prep=False,
                          run_seas_phyto=False,
                          run_match_phyto1=False,
                          run_eval_phyto1=False,
                          run_match_phyto2=False,
                          run_eval_phyto2=False,
                          run_eval_bt1=True,
                          run_nutr=False,
                          run_zp_eval=False):
    if run_prep:
        print("=== Running Data Prep Step ===")
        data_prep.run_ecospace_prep()

    if run_seas_phyto:
        print("=== Running Seas Phyto (McEwan) ===")
        phyto_seas.run_eval_2d_phyto()

    if run_match_phyto1:
        print("=== Running Phyto Match 1 (Nemcek) ===")
        phyto_match_1.run_eval_2a_nemcek_match()

    if run_eval_phyto1:
        print("=== Running Phyto Eval 1 (Nemcek)"
              " ===")
        phyto_eval_1.run_eval_2b_nemcek_eval()

    if run_match_phyto2:
        print("=== Running Phyto Match 2 (QU39) ===")
        phyto_match_2.run_qu39_match()

    if run_eval_phyto2:
        print("=== Running Phyto Eval 2 (QU39) ===")
        phyto_eval_2.run_qu39_eval()

    if run_eval_bt1:
        print("=== Running BT Eval 1 (CSOG) ===")
        bt_eval_1.run_bt_eval()

    if run_nutr:
        print("=== Running nutrient eval ===")
        # nu_eval.

    if run_zp_eval:
        print("=== Running zoop eval ===")
        zp_eval.run_zoop_eval()

if __name__ == "__main__":
    run_ecospace_pipeline()

