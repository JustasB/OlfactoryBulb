
if __name__ == '__main__':

    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params
    print 'New folder:', params['folder_name']
    ok = raw_input('\nContinue to create this folder structure? Parameters therein will be overwritten\n\ty / Y / blank = OK; anything else --> exit\n')
    if not ((ok == '') or (ok.capitalize() == 'Y')):
        print 'quit'
        exit(1)

    input_folders = ['PaperData/Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/', 
            'Cluster_ORnoise0.05_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/', 
            'Cluster_ORnoise0.10_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']
#            'Cluster_ORnoise0.15_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']
#            'Cluster_ORnoise0.2_nGlom40_nHC12_nMC30_vqOvrlp4_np50_fullSystem/']
#            'Cluster_ORnoise0.25_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']
#            'Cluster_ORnoise0.30_nGlom40_nHC12_nMC30_vqOvrlp4_np50_FullSystem/']

    n_pattern_sets = len(inupt_folders)
