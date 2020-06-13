# MOCKS for autodoc
import quantities as pq
if pq.mV.__class__.__module__ == 'sphinx.ext.autodoc.importer':
    pq.mV = pq.ms = pq.Hz = pq.nA = 1.0
# END MOCKS

from abc import abstractmethod

from neuronunit import capabilities as ncap
from neuronunit.tests.base import VmTest
from olfactorybulb.neuronunit.tests.utilities import get_APs, cache
from sciunit import capabilities as scap

from olfactorybulb.neuronunit import capabilities as obncap
from olfactorybulb.neuronunit.tests import publications

SHOW_ERRORS = False

class OlfactoryBulbCellTest(VmTest):

    @abstractmethod
    def generate_prediction_nocache(self, model):
        pass

    def generate_prediction(self, model):
        # import pydevd
        # pydevd.settrace('192.168.0.100', port=4200, suspend=False)

        result = self.fetch_cached(model)

        if result is None:

            # Check that self has all the required properties
            self.check_required_properties()

            # Perform the uncached test
            try:
                result = self.generate_prediction_nocache(model)
            except:
                import traceback
                result = traceback.format_exc()

                if SHOW_ERRORS:
                    print(result)

            # Store result in cache
            self.store_in_cache(model, result)

        return result

    def check_required_properties(self):
        if hasattr(self, "required_properties"):
            for prop in self.required_properties:
                if not hasattr(self, prop):
                    raise Exception("Property '" + prop + "' not found. Make sure the property is declared either in the"
                                                         " generic test class or in the publication class.")

    def fetch_cached(self, model):
        return cache.get(self.get_hash(model))

    def store_in_cache(self, model, result):
        cache.store(self.get_hash(model), result)

    def get_hash(self, model):
        # The cache key is a hash of the model and the test - we want to store the model-test_result combination
        model_hash = model.__hash__()
        self_hash = self.__hash__()
        return hash((model_hash, self_hash))

    def __hash__(self):
        return hash(self.__class__.__name__)

    def get_dependent_prediction(self, dependent_test_class_generic, model):

        # import pydevd
        # pydevd.settrace('192.168.0.100', port=4200)

        mro = self.__class__.mro()

        if len(mro) < 4:
            raise Exception("The test should be a class that inherits from an publications class"
                            "AND from a generic tests class, in that order. E.g. "
                            "'class MyTest(UrbanBurton2014, InputResistanceTest):'")

        # Create a temp class that inherits from the generic test and from the specific publication
        # Aways first parent class (by convention and to preserve inheritance)
        publication_class = mro[1]


        if not issubclass(publication_class, publications.BasePublication):
            raise Exception("The first parent class '"+str(publication_class)+"' of the test should be a publication class. E.g. 'class MyTest(UrbanBurton2014, InputResistanceTest):'")


        if not issubclass(dependent_test_class_generic, OlfactoryBulbCellTest):
            raise Exception("The second parent class '"+dependent_test_class_generic.__class__.__name__+"' of the test should be a class that inherits from OlfactoryBulbCellTest. E.g. 'class MyTest(UrbanBurton2014, InputResistanceTest):'")


        # Use SomeTestSomeAuthor1984 class name form - as descibed in BasePublication
        dependent_test_class_name = dependent_test_class_generic.__name__ + publication_class.__name__

        # Create the type dynamically
        dependent_test_class = type(dependent_test_class_name,
                                    (publication_class, dependent_test_class_generic),
                                    {})

        # Instantiate the dynamic class
        dependent_test = dependent_test_class()

        # Get the prediction (from cache if there)
        return dependent_test.generate_prediction(model)



class OlfactoryBulbCellSpikeTest(OlfactoryBulbCellTest):

    required_capabilities = (ncap.ReceivesSquareCurrent,
                             ncap.ProducesMembranePotential,
                             scap.Runnable,
                             obncap.SupportsSettingTemperature,
                             obncap.SupportsSettingStopTime)

    def get_aps(self, voltage):
        return get_APs(voltage, self.ss_delay, self.threshold_method)