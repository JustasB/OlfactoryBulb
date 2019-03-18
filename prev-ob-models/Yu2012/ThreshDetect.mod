NEURON {
     POINT_PROCESS ThreshDetect

     RANGE vthresh
 }

UNITS {

	(mV) = (millivolt)

}
 PARAMETER {
        vthresh= -40 (mV) 
         
}
 ASSIGNED { v (mV) }
 INITIAL {
     net_send(0, 1)
 }
 NET_RECEIVE(w) {
     if (flag == 1) {
        :printf ("flag =1")
        WATCH (v > vthresh) 2
     }else if (flag == 2) {
        :printf ("flag =2")
        net_event(t)
        
        :printf ("v=%g",v)

     }
 }
