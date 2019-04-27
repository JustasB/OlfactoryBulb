: Adds the function flushf() to hoc. This is used for the
: cell-creation progress meter.

NEURON {
	SUFFIX nothing
}

PROCEDURE flushf() {
  VERBATIM
    fflush(NULL);
  ENDVERBATIM
}
