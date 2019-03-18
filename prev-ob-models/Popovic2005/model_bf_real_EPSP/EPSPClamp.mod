TITLE svclmp.mod

COMMENT

Single electrode Voltage clamp with three levels
------------------------------------------------

Series Resistance added; backards compatible, except parameters 
e0,vo0,vi0,gain,rstim,tau1,tau2 that no longer exist

Clamp is on at time 0, and off at time dur[0]+dur[1]+dur[2]. When clamp is off
the injected current is 0.  The clamp levels are amp[0], amp[1], amp[2].  i is
the injected current, vc measures the control voltage) Do not insert several
instances of this model at the same location in order to make level changes.
That is equivalent to independent clamps and they will have incompatible
internal state values.

The electrical circuit for the clamp is exceedingly simple:

        rs           Rin
vc ---'\/\/`---o---'\/\/`---o
               |            |
               |____| |_____|
                    | |
                     Cm

Note that since this is an electrode current model v refers to the internal
potential which is equivalent to the membrane potential v when there is no
extracellular membrane mechanism present but is v+vext when one is present. 
Also since i is an electrode current, positive values of i depolarize the
cell. (Normally, positive membrane currents are outward and thus hyperpolarize
the cell)

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 200

NEURON {
        POINT_PROCESS EPSPClamp
        ELECTRODE_CURRENT i
        RANGE delay, rs, vc, i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (micromho)
}


PARAMETER {
        v (mV)
        rs = 1 (megohm)		: series resistance
}

ASSIGNED {
        i (nA)
        vc (mV)
        ic (nA)
	  delay (ms)
        on
}

INITIAL {
        on = 0
}

BREAKPOINT {
        SOLVE vstim
        if (on) {
                i = (vc - v)/rs
        }else{
                i = 0
        }
}

PROCEDURE vstim() {
        on = 1
        if (t < delay) { vc = -55
        }else if (t < delay+0.0) { vc = -55.0000
        }else if (t < delay+0.5) { vc = -54.9824
        }else if (t < delay+1.0) { vc = -53.5728
        }else if (t < delay+1.5) { vc = -50.0865
        }else if (t < delay+2.0) { vc = -46.7973
        }else if (t < delay+2.5) { vc = -44.5711
        }else if (t < delay+3.0) { vc = -43.3592
        }else if (t < delay+3.5) { vc = -42.6923
        }else if (t < delay+4.0) { vc = -42.1744
        }else if (t < delay+4.5) { vc = -41.9764
        }else if (t < delay+5.0) { vc = -42.0832
        }else if (t < delay+5.5) { vc = -42.1900
        }else if (t < delay+6.0) { vc = -42.6444
        }else if (t < delay+6.5) { vc = -43.0820
        }else if (t < delay+7.0) { vc = -43.5035
        }else if (t < delay+7.5) { vc = -43.9092
        }else if (t < delay+8.0) { vc = -44.2995
        }else if (t < delay+8.5) { vc = -44.6750
        }else if (t < delay+9.0) { vc = -45.0360
        }else if (t < delay+9.5) { vc = -45.3829
        }else if (t < delay+10.0) { vc = -45.7162
        }else if (t < delay+10.5) { vc = -46.0363
        }else if (t < delay+11.0) { vc = -46.3435
        }else if (t < delay+11.5) { vc = -46.6382
        }else if (t < delay+12.0) { vc = -46.9209
        }else if (t < delay+12.5) { vc = -47.1918
        }else if (t < delay+13.0) { vc = -47.4514
        }else if (t < delay+13.5) { vc = -47.7000
        }else if (t < delay+14.0) { vc = -47.9380
        }else if (t < delay+14.5) { vc = -48.1657
        }else if (t < delay+15.0) { vc = -48.3834
        }else if (t < delay+15.5) { vc = -48.5915
        }else if (t < delay+16.0) { vc = -48.7902
        }else if (t < delay+16.5) { vc = -48.9800
        }else if (t < delay+17.0) { vc = -49.1610
        }else if (t < delay+17.5) { vc = -49.3337
        }else if (t < delay+18.0) { vc = -49.4982
        }else if (t < delay+18.5) { vc = -49.6549
        }else if (t < delay+19.0) { vc = -49.8041
        }else if (t < delay+19.5) { vc = -49.9460
        }else if (t < delay+20.0) { vc = -50.0810
        }else if (t < delay+20.5) { vc = -50.2091
        }else if (t < delay+21.0) { vc = -50.3308
        }else if (t < delay+21.5) { vc = -50.4463
        }else if (t < delay+22.0) { vc = -50.5558
        }else if (t < delay+22.5) { vc = -50.6595
        }else if (t < delay+23.0) { vc = -50.7577
        }else if (t < delay+23.5) { vc = -50.8506
        }else if (t < delay+24.0) { vc = -50.9384
        }else if (t < delay+24.5) { vc = -51.0213
        }else if (t < delay+25.0) { vc = -51.0996
        }else if (t < delay+25.5) { vc = -51.1733
        }else if (t < delay+26.0) { vc = -51.2429
        }else if (t < delay+26.5) { vc = -51.3083
        }else if (t < delay+27.0) { vc = -51.3699
        }else if (t < delay+27.5) { vc = -51.4277
        }else if (t < delay+28.0) { vc = -51.4821
        }else if (t < delay+28.5) { vc = -51.5330
        }else if (t < delay+29.0) { vc = -51.5807
        }else if (t < delay+29.5) { vc = -51.6254
        }else if (t < delay+30.0) { vc = -51.6672
        }else if (t < delay+30.5) { vc = -51.7063
        }else if (t < delay+31.0) { vc = -51.7428
        }else if (t < delay+31.5) { vc = -51.7768
        }else if (t < delay+32.0) { vc = -51.8085
        }else if (t < delay+32.5) { vc = -51.8380
        }else if (t < delay+33.0) { vc = -51.8655
        }else if (t < delay+33.5) { vc = -51.8910
        }else if (t < delay+34.0) { vc = -51.9147
        }else if (t < delay+34.5) { vc = -51.9367
        }else if (t < delay+35.0) { vc = -51.9571
        }else if (t < delay+35.5) { vc = -51.9760
        }else if (t < delay+36.0) { vc = -51.9935
        }else if (t < delay+36.5) { vc = -52.0097
        }else if (t < delay+37.0) { vc = -52.0247
        }else if (t < delay+37.5) { vc = -52.0386
        }else if (t < delay+38.0) { vc = -52.0515
        }else if (t < delay+38.5) { vc = -52.0634
        }else if (t < delay+39.0) { vc = -52.0745
        }else if (t < delay+39.5) { vc = -52.0848
        }else if (t < delay+40.0) { vc = -52.0944
        }else if (t < delay+40.5) { vc = -52.1033
        }else if (t < delay+41.0) { vc = -52.1116
        }else if (t < delay+41.5) { vc = -52.1194
        }else if (t < delay+42.0) { vc = -52.1268
        }else if (t < delay+42.5) { vc = -52.1337
        }else if (t < delay+43.0) { vc = -52.1403
        }else if (t < delay+43.5) { vc = -52.1466
        }else if (t < delay+44.0) { vc = -52.1526
        }else if (t < delay+44.5) { vc = -52.1584
        }else if (t < delay+45.0) { vc = -52.1640
        }else if (t < delay+45.5) { vc = -52.1694
        }else if (t < delay+46.0) { vc = -52.1748
        }else if (t < delay+46.5) { vc = -52.1801
        }else if (t < delay+47.0) { vc = -52.1853
        }else if (t < delay+47.5) { vc = -52.1905
        }else if (t < delay+48.0) { vc = -52.1957
        }else if (t < delay+48.5) { vc = -52.2010
        }else if (t < delay+49.0) { vc = -52.2063
        }else if (t < delay+49.5) { vc = -52.2117
        }else if (t < delay+50.0) { vc = -52.2171
        }else if (t < delay+50.5) { vc = -52.2227
        }else if (t < delay+51.0) { vc = -52.2284
        }else if (t < delay+51.5) { vc = -52.2341
        }else if (t < delay+52.0) { vc = -52.2401
        }else if (t < delay+52.5) { vc = -52.2461
        }else if (t < delay+53.0) { vc = -52.2523
        }else if (t < delay+53.5) { vc = -52.2586
        }else if (t < delay+54.0) { vc = -52.2651
        }else if (t < delay+54.5) { vc = -52.2717
        }else if (t < delay+55.0) { vc = -52.2785
        }else if (t < delay+55.5) { vc = -52.2853
        }else if (t < delay+56.0) { vc = -52.2924
        }else if (t < delay+56.5) { vc = -52.2995
        }else if (t < delay+57.0) { vc = -52.3068
        }else if (t < delay+57.5) { vc = -52.3142
        }else if (t < delay+58.0) { vc = -52.3217
        }else if (t < delay+58.5) { vc = -52.3292
        }else if (t < delay+59.0) { vc = -52.3369
        }else if (t < delay+59.5) { vc = -52.3446
        }else if (t < delay+60.0) { vc = -52.3524
        }else if (t < delay+60.5) { vc = -52.3602
        }else if (t < delay+61.0) { vc = -52.3681
        }else if (t < delay+61.5) { vc = -52.3759
        }else if (t < delay+62.0) { vc = -52.3838
        }else if (t < delay+62.5) { vc = -52.3916
        }else if (t < delay+63.0) { vc = -52.3994
        }else if (t < delay+63.5) { vc = -52.4072
        }else if (t < delay+64.0) { vc = -52.4148
        }else if (t < delay+64.5) { vc = -52.4224
        }else if (t < delay+65.0) { vc = -52.4298
        }else if (t < delay+65.5) { vc = -52.4372
        }else if (t < delay+66.0) { vc = -52.4444
        }else if (t < delay+66.5) { vc = -52.4514
        }else if (t < delay+67.0) { vc = -52.4582
        }else if (t < delay+67.5) { vc = -52.4649
        }else if (t < delay+68.0) { vc = -52.4713
        }else if (t < delay+68.5) { vc = -52.4774
        }else if (t < delay+69.0) { vc = -52.4833
        }else if (t < delay+69.5) { vc = -52.4890
        }else if (t < delay+70.0) { vc = -52.4943
        }else if (t < delay+70.5) { vc = -52.4994
        }else if (t < delay+71.0) { vc = -52.5041
        }else if (t < delay+71.5) { vc = -52.5084
        }else if (t < delay+72.0) { vc = -52.5124
        }else if (t < delay+72.5) { vc = -52.5161
        }else if (t < delay+73.0) { vc = -52.5193
        }else if (t < delay+73.5) { vc = -52.5222
        }else if (t < delay+74.0) { vc = -52.5246
        }else if (t < delay+74.5) { vc = -52.5266
        }else if (t < delay+75.0) { vc = -52.5282
        }else if (t < delay+75.5) { vc = -52.5293
        }else if (t < delay+76.0) { vc = -52.5300
        }else if (t < delay+76.5) { vc = -52.5302
        }else if (t < delay+77.0) { vc = -52.5299
        }else if (t < delay+77.5) { vc = -52.5292
        }else if (t < delay+78.0) { vc = -52.5280
        }else if (t < delay+78.5) { vc = -52.5264
        }else if (t < delay+79.0) { vc = -52.5242
        }else if (t < delay+79.5) { vc = -52.5216
        }else if (t < delay+80.0) { vc = -52.5184
        }else if (t < delay+80.5) { vc = -52.5149
        }else if (t < delay+81.0) { vc = -52.5108
        }else if (t < delay+81.5) { vc = -52.5063
        }else if (t < delay+82.0) { vc = -52.5013
        }else if (t < delay+82.5) { vc = -52.4959
        }else if (t < delay+83.0) { vc = -52.4901
        }else if (t < delay+83.5) { vc = -52.4838
        }else if (t < delay+84.0) { vc = -52.4771
        }else if (t < delay+84.5) { vc = -52.4701
        }else if (t < delay+85.0) { vc = -52.4627
        }else if (t < delay+85.5) { vc = -52.4549
        }else if (t < delay+86.0) { vc = -52.4469
        }else if (t < delay+86.5) { vc = -52.4385
        }else if (t < delay+87.0) { vc = -52.4299
        }else if (t < delay+87.5) { vc = -52.4210
        }else if (t < delay+88.0) { vc = -52.4119
        }else if (t < delay+88.5) { vc = -52.4027
        }else if (t < delay+89.0) { vc = -52.3934
        }else if (t < delay+89.5) { vc = -52.3840
        }else if (t < delay+90.0) { vc = -52.3745
        }else if (t < delay+90.5) { vc = -52.3650
        }else if (t < delay+91.0) { vc = -52.3556
        }else if (t < delay+91.5) { vc = -52.3463
        }else if (t < delay+92.0) { vc = -52.3371
        }else if (t < delay+92.5) { vc = -52.3282
        }else if (t < delay+93.0) { vc = -52.3195
        }else if (t < delay+93.5) { vc = -52.3112
        }else if (t < delay+94.0) { vc = -52.3033
        }else if (t < delay+94.5) { vc = -52.2958
        }else if (t < delay+95.0) { vc = -52.2889
        }else if (t < delay+95.5) { vc = -52.2827
        }else if (t < delay+96.0) { vc = -52.2771
        }else if (t < delay+96.5) { vc = -52.2722
        }else if (t < delay+97.0) { vc = -52.2683
        }else if (t < delay+97.5) { vc = -52.2653
        }else if (t < delay+98.0) { vc = -52.2633
        }else if (t < delay+98.5) { vc = -52.2625
        }else if (t < delay+99.0) { vc = -52.2629
        }else if (t < delay+99.5) { vc = -52.2647
        }else {
                vc = -52.2647
                on = 0
        }
        if (on) {
        }else{
                ic = 0
        }
        VERBATIM
        return 0;
        ENDVERBATIM
}
