# MOCKS for autodoc
import quantities as pq
if pq.mV.__class__.__module__ == 'sphinx.ext.autodoc.importer':
    pq.mV = pq.ms = pq.Hz = pq.nA = pq.pA = 1
# END MOCKS