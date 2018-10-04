from database import *
import numpy as np


def pool_experimental_measurements(property_id):
    """
    Aggregates the distributions of measurments of the given property and saves
    the resulting aggregate mean, std, and n to property table

    :param property_id: the id of the property to aggregate
    """

    measurements = Measurement\
        .select()\
        .where(Measurement.property == property_id)

    simulated_values = []

    for m in measurements:
        simulated_values += list(np.random.normal(m.mean, m.std, m.n))

    prop = Property\
        .get(Property.id == property_id)

    prop.mean = np.mean(simulated_values)
    prop.std = np.std(simulated_values)
    prop.n = len(simulated_values)
    prop.save()

    print("Pooled %d measurements of %s with mean: %s std: %s n: %d" %
        (len(measurements), prop.id, prop.mean, prop.std, prop.n))


def pick_from_pooled_experimental_distribution(property_id, type='float', save=True):
    """
    Picks one value, at random, from the property's aggregate distribution

    :param property_id: the id of the property, whose distribution will be used to pick
    :param type: the type of value to return e.g. use 'int' to round final result to integer
    :param save: whether to save the picked result to the database
    :return:
    """
    prop = Property \
        .get(Property.id == property_id)

    picked = np.random.normal(prop.mean, prop.std)

    if type == 'int':
        picked = int(round(picked))

    prop.picked_value = picked

    if save:
        prop.save()

    return picked



