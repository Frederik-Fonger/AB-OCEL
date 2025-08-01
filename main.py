import sampling.manager
import evaluation.OCEL_evaluation



if __name__ == '__main__':

    # -------------- Sampling -----------------

    sampling = sampling.manager.OCEL_sampling()
    sampling.apply('data/02_p2p.xml', 0.05, sampling_algo="AB", file_type="XML", create_folder_with_sample_and_meta_data = True )
    


    # --------------  Evaluation  --------------
    # eval = evaluation.OCEL_evaluation.Evaluation()
    
    # ------ fitness -----
    # eval.eval_from_file('data/lrms collection 02 p2p/AB_OCEL/AB_0_3/02_p2p_sample_AB_0.3.xmlocel' ,'data/lrms collection 02 p2p/AB_OCEL/AB_0_3/02_p2p_sample_AB_0.3.xmlocel', token=False, output_folder='/output/eval')

    # --------- MAE and coverage --------
    # eval.calculate_MAE_and_coverage(sample='data/lrms collection 02 p2p/AB_OCEL/AB_0_3/02_p2p_sample_AB_0.3.xmlocel', original_log='data/02_p2p.xml')
