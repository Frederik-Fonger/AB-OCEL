import time

import pm4py
import pandas as pd
import sampling.manager as smpl
from ocpa.objects.log.importer.ocel2.xml import factory as ocel_import_factory
from ocpa.algo.discovery.ocpn import algorithm as ocpn_discovery_factory
from ocpa.algo.conformance.precision_and_fitness import evaluator as quality_measure_factory
# from ocpa.algo.enhancement.token_replay_based_performance import algorithm as token_based_replay_algorithm
from ocpa_main.ocpa.algo.conformance.token_based_replay import algorithm as token_based_replay_algorithm
class Evaluation:


    def __init__(self):
        """
        Stores the variable need while evaluating. The variable are set, as soon as there are available.
        """
        self.sample_pm4py = None
        self.sample_ocpa = None
        self.ocel_pm4py = None
        self.ocel_pm4py_orginal = None
        self.ocel_ocpa = None
        self.ocel_ocpa_orginal = None
        self.sampling_algo = None
        self.sampling_ratio = None
        self.quality_metrics = None
        self.token_replay = False


    def import_ocel(self, path):
        self.ocel_pm4py = pm4py.read.read_ocel2_xml(path)
        self.ocel_pm4py_orginal = pm4py.read.read_ocel2_xml(path)
        self.ocel_ocpa = ocel_import_factory.apply(path)
        self.ocel_ocpa_orginal = ocel_import_factory.apply(path)


    def update_sample_from_pm4py_to_ocpa(self):
        """
        This function updates changes from the pm4py-ocel object to the ocpa-ocel object. This is done by save and load of a temp XML file.
        :return:
        """
        pm4py.write.write_ocel2_xml(self.sample_pm4py, "temp.xmlocel")
        self.sample_ocpa = ocel_import_factory.apply("temp.xmlocel")

    def update_from_ocpa_to_pm4py(self):
        # Todo this function
        pass

    def calculate_quality_metrics(self, ocel, ocpn, output_folder):

        if self.token_replay:
            metrics_dict = token_based_replay_algorithm.apply(ocel, ocpn)
        else:
            precision, fitness = quality_measure_factory.apply(ocel, ocpn, output_folder)
            metrics_dict = {
                'fitness' : fitness,
                'precision': precision
            }

        self.quality_metrics = metrics_dict


    def sample_IM_metrics(self, path, sample_ratio, sampling_algo="AB", connectivity_threshold=0.5):

        # import log
        Evaluation.import_ocel(self,path)

        # my sampling
        sampling = smpl.OCEL_sampling()
        self.sample_pm4py = sampling.striped_sampling(self.ocel_pm4py, self.ocel_ocpa, sample_ratio, sampling_algo, connectivity_threshold)
        Evaluation.update_sample_from_pm4py_to_ocpa(self)

        #IM
        # ocel = ocel_import_factory.apply('result.xmlocel')
        print("reults-xml imported")
        ocpn = ocpn_discovery_factory.apply(self.sample_ocpa, parameters={"debug": True})
        print("Petri net mined")
        # metrics on model and sample

        Evaluation.calculate_quality_metrics(self, self.sample_ocpa, ocpn)

        print(self.quality_metrics)

    def eval_from_file(self, path_log_for_model, path_log_for_eval, token, output_folder):

        print("Log for model: " + path_log_for_model)
        print("Log for eval:  " + path_log_for_eval)
        print("Output folder: " + output_folder)

        #IM
        self.token_replay = token
        print("Path for ocel factory: " + path_log_for_model) 
        self.sample_ocpa = ocel_import_factory.apply(path_log_for_model)
        print("reults-xml imported")
        ocpn = ocpn_discovery_factory.apply(self.sample_ocpa, parameters={"debug": True})
        print("Petri net mined")
        # metrics on model and sample

        self.ocel_ocpa_orginal = ocel_import_factory.apply(path_log_for_eval)
        start_time = time.time()
        Evaluation.calculate_quality_metrics(self, self.ocel_ocpa_orginal, ocpn, output_folder)
        end_time = time.time()
        print("Time for eval was: " + str((end_time - start_time)) + "seconds ")
        print(self.quality_metrics)

        data = {'Info': ['path_log_for_model','path_log_for_eval','token','metrics' ],
                'Value': [path_log_for_model, path_log_for_eval, token, self.quality_metrics]}
        df = pd.DataFrame(data)
        df.to_csv(output_folder + "output.csv", index=False)