import params

def mgid2glom(gid):
  return int((gid-params.gid_mitral_begin)/params.Nmitral_per_glom)
