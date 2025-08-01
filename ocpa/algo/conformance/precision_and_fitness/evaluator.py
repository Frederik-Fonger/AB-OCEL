from ocpa.algo.conformance.precision_and_fitness.variants import replay_context
import ocpa.algo.conformance.precision_and_fitness.utils as utils
import time

def apply(ocel,ocpn,output_folder, contexts=None,bindings=None):
    '''
    Calculation precision and fitness for an object-centric Petri net with respect to an object-centric event log. The
    measures are calculated according to replaying the event log and checking enabled and executed behavior. Contexts and
    bindings can be pre-computed and passed to the method to save computation time upon multiple calling. If not given,
    contexts and binding wil be calculated.

    :param ocel: Object-centric event log
    :type ocel: :class:`OCEL <ocpa.objects.log.ocel.OCEL>`

    :param ocpn: Object-centric Petri net
    :type ocpn: :class:`OCPN <ocpa.objects.oc_petri_net.obj.ObjectCentricPetriNet>`

    :param contexts: multiset of previously executed traces of activities for each event (can be computed by calling :func:`the corresponding function <ocpa.algo.evaluation.precision_and_fitness.utils.calculate_contexts_and_bindings>`)
    :type contexts: Dict

    :param bindings: bindings for each event (can be computed by calling :func:`the corresponding function <ocpa.algo.evaluation.precision_and_fitness.utils.calculate_contexts_and_bindings>`)
    :type bindings: Dict

    :return: precision, fitness
    :rtype: float, float

    '''
    object_types = ocel.object_types
    print("Calculating context")
    start_time = time.time()
    if contexts == None or bindings == None:
        contexts, bindings = (utils.calculate_contexts_and_bindings(ocel))
    end_time = time.time()
    print("Time for Calculating context was: " + str((end_time - start_time)) + "seconds ")
    print("calcualting en_l")
    start_time = time.time()
    en_l =  replay_context.enabled_log_activities(ocel.log, contexts)
    end_time = time.time()
    print("Time for calcualting en_l was: " + str((end_time - start_time)) + "seconds ")
    print("calcualting en_m")
    start_time = time.time()
    en_m =  replay_context.enabled_model_activities_multiprocessing(contexts, bindings, ocpn, object_types, output_folder)
    end_time = time.time()
    print("Time for calcualting en_m was: " + str((end_time - start_time)) + "seconds ")
    print("Calculating precision and fitness")
    start_time = time.time()
    precision, skipped_events, fitness =  replay_context.calculate_precision_and_fitness(ocel.log, contexts, en_l, en_m)
    end_time = time.time()
    print("Time for Calculating precision and fitness was: " + str((end_time - start_time)) + "seconds ")
    return precision, fitness
    