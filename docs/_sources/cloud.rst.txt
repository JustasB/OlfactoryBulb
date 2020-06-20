Running Simulations in the Cloud (AWS)
======================================

To accelerate simulations, the model can be executed on cloud virtual machines. For example, as of this writing,
a 96 Core machine (``c5.24xlarge``) can be rented from Amazon Web Services for ~$1.5/hr (spot price). The basic usage cycle involves the
following steps:

 - Request a 'spot' instance from AWS (spot instances are ~70% less expensive)
 - SSH into the instance, download the model docker image, run the simulation
 - Upload simulation results for later analysis
 - Shut down instance (per-second billing)

A similar process can be followed with other cloud providers (GCP, Azure, etc...).

==========================
Private Cloud MPI Clusters
==========================

Similarly, multiple instances could be started to form a private MPI cluster, and the simulation could be distributed across multiple instances. See `AWS ParallelCluster <https://docs.aws.amazon.com/parallelcluster/latest/ug/what-is-aws-parallelcluster.html>`_.



