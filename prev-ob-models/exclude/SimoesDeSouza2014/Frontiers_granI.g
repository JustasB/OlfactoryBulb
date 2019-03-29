// genesis 2.3

include defaults

int hflag = 1    // use hsolve if hflag = 1
int chanmode = 1

float tmax = 0.1               // simulation time in sec
float dt = 1e-7              // simulation time step in sec
setclock  0  {dt}               // set the simulation clock
float injcurr = 0

setclock 1 1e-4 //clock for plots


// Create a library of prototype elements to be used by the cell reader


/* file for standard compartments */
include compartments 
//* file for Hodgkin-Huxley Squid Na and K channels */
include hhchan 

/* file for Upi's mitral cell channels */
include mitchan 
include newbulbchan

/* file for Upi's mitral cell synaptic channels */
include mitsynC2 // for now use channelC2 version
// file for Fabio's prototype
include fabioproto

/************************************************************************
**  2	Invoking functions to make prototypes in the /library element
************************************************************************/


//create neutral /library

disable /library

	pushe /library


	/* Make the standard types of compartments  */

	make_cylind_compartment		/* makes "compartment" */
	make_sphere_compartment		/* makes "compartment_sphere" */
	make_cylind_symcompartment	/* makes "symcompartment" */
	make_sphere_symcompartment	/* makes "symcompartment_sphere" */

	/* These are some standard channels used in .p files */
	make_Na_squid_hh		/* makes "Na_squid_hh" */
	make_K_squid_hh			/* makes "K_squid_hh" */
	make_Na_mit_hh			/* makes "Na_mit_hh" */
	make_K_mit_hh			/* makes "K_mit_hh" */

make_LCa3_mit_usb
make_Na_rat_smsnn	// Na current
make_Na2_rat_smsnn
make_KA_bsg_yka
make_KM_bsg_upi
make_K_mit_usb  // K-current     
make_K2_mit_usb
make_K_slow_usb
make_Na_mit_usb
make_Na2_mit_usb
make_glu_mit_usb
make_GABA_1_mit_usb
make_GABA_2_mit_usb
make_glu_gran_usb
make_glu_pg_usb
make_olf_receptor
make_spike
make_Kca_mit_usb
make_Ca_mit_conc

make_GABA_A //fabioproto
make_AMPA_NMDA // fabioproto
make_Ca_conc /fabioproto



	/* These are some synaptic channels for the mitral cell */
	make_glu_mit_upi		/* makes "glu_mit_upi" */
	make_GABA_mit_upi		/* makes "GABA_mit_upi" */

	pope

enable /library

//===============================
//      Function Definitions
//===============================

function step_tmax
    step {tmax} -time
end

//===============================
//    Graphics Functions
//===============================


function make_control
    create xform /control [10,50,300,200]
    create xlabel /control/label -hgeom 50 -bg cyan -label "CONTROL PANEL"
    create xbutton /control/RESET -wgeom 25%       -script reset
    create xbutton /control/RUN  -xgeom 0:RESET -ygeom 0:label -wgeom 19% \
         -script step_tmax
    create xbutton /control/STOP  -xgeom 0:RUN -ygeom 0:label -wgeom 19% \
         -script stop
    create xbutton /control/QUIT -xgeom 0:STOP -ygeom 0:label -wgeom 19% \
        -script quit
    create xbutton /control/STIM -xgeom 0:QUIT -ygeom 0:label -wgeom 19% \
        -script stim
    create xdialog /control/Injection -label "Injection (amperes)" \
		-value {injcurr}  -script "set_inject <widget>"
    create xdialog /control/stepsize -title "dt (sec)" -value {dt} \
                -script "change_stepsize <widget>"
    create xdialog /control/tempo -label "Tmax (sec)" \
                -value {tmax} -script "set_time <widget>"
    create xtoggle /control/overlay  \
           -script "overlaytoggle <widget>"
    setfield /control/overlay offlabel "Overlay OFF" onlabel "Overlay ON" state 0
    xshow /control
end

function make_Vmgraph
    float vmin = -0.100
    float vmax = 0.05
    create xform /data [265,50,350,350]
    create xlabel /data/label -hgeom 10% -label "Granule cell "
    create xgraph /data/voltage  -hgeom 90%  -title "Membrane Potential"
    setfield ^ XUnits sec YUnits Volts
    setfield ^ xmax {tmax} ymin {vmin} ymax {vmax}
    xshow /data
    useclock /data/## 1 
end

function make_Gkgraph
    float gmin = 0
    float gmax = 1e-9
    create xform /data2 [265,50,350,350]
    create xlabel /data2/label -hgeom 10% -label "Granule cell"
    create xgraph /data2/conductance  -hgeom 90%  -title "Conductances"
    setfield ^ XUnits sec YUnits Siemens
    setfield ^ xmax {tmax} ymin {gmin} ymax {gmax}
    xshow /data2
    useclock /data2/## 1 
end

function make_xcell
    create xform /cellform [620,50,400,400]
    create xdraw /cellform/draw [0,0,100%,100%]
    setfield /cellform/draw xmin -3e-3 xmax 3e-3 ymin -3e-3 ymax 3e-3 zmin -3e-2 zmax 3e-2 transform p
    xshow /cellform
    echo creating xcell
    create xcell /cellform/draw/cell
    setfield /cellform/draw/cell colmin -0.1 colmax 0.1 \
        path /cell/##[TYPE=compartment] field Vm \
        script "echo <w> <v>"
    useclock /cellform/## 1 
end


function set_time(dialog)
    str dialog
    tmax={getfield {dialog} value}
end

function set_inject(dialog)
    str dialog
    setfield /cell/soma inject {getfield {dialog} value}
end

function change_stepsize(dialog)
   str dialog
   dt =  {getfield {dialog} value}
   setclock 0 {dt}
   echo "Changing step size to "{dt}
end

function stim
step 0.010 -time
echo "Firing synapse"
setfield presyn z {1/{getclock 0}}
step 1
setfield presyn z 0
step 1
step 0.9 -time
echo "End of simulation"
end

// Use of the wildcard sets overlay field for all graphs
function overlaytoggle(widget)
    str widget
    setfield /##[TYPE=xgraph] overlay {getfield {widget} state}
end

//===============================
//         Main Script
//===============================
// Build the cell from a parameter file using the cell reader

readcell granI_fabio.p /cell -hsolve 

setfield /cell path "/cell/##[][TYPE=compartment]"
setfield /cell chanmode {chanmode}
call /cell SETUP
setmethod 11
echo "Using hsolve"
reset

setfield /cell/soma inject {injcurr}


//setting up ca_fraction and Enmda
function set_Cafraction
//fractional Ca2+ (Schneggenburger et al., 1996)
float Caout=2 //mM extracelular calcium
float M=155  //mM extracelular monovalent ion concentration
float PCaPM=3.6  //permeability ration of Ca2+ over monovalent ions
float T=298.15 //273.15+25 Celcius Temperature
float R=8.314 //J.K-1.mol-1 gas constant
float F=96485 //C.mol-1 Farady constant
float V
float f

//GHK approach for Enmda
float Enmda
float a
a={{4*M*{M+{4*PCaPM*Caout}}**0.5}/{2*M}}
Enmda= {{R*T}/F}*{log {a}}
setfield /cell/periph12[17]/head[17]/AMPA_NMDA Ek Enmda
setfield /cell/periph11[20]/head[20]/AMPA_NMDA Ek Enmda
setfield /cell/periph21[5]/head[5]/AMPA_NMDA Ek Enmda
setfield /cell/periph22[8]/head[8]/AMPA_NMDA Ek Enmda

V={getfield /cell/periph12[17]/head[17] Vm}
f={Caout}/{Caout+{1/PCaPM}*{M/4}*{1-{exp {2*V*F/{R*T}}}}}
setfield /cell/periph12[17]/head[17]/AMPA_NMDA/Ca_fraction x {f}

V={getfield /cell/periph11[20]/head[20] Vm}
f={Caout}/{Caout+{1/PCaPM}*{M/4}*{1-{exp {2*V*F/{R*T}}}}}
setfield /cell/periph11[20]/head[20]/AMPA_NMDA/Ca_fraction x {f}

V={getfield /cell/periph21[5]/head[5] Vm}
f={Caout}/{Caout+{1/PCaPM}*{M/4}*{1-{exp {2*V*F/{R*T}}}}}
setfield /cell/periph21[5]/head[5]/AMPA_NMDA/Ca_fraction x {f}

V={getfield /cell/periph22[8]/head[8] Vm}
f={Caout}/{Caout+{1/PCaPM}*{M/4}*{1-{exp {2*V*F/{R*T}}}}}
setfield /cell/periph22[8]/head[8]/AMPA_NMDA/Ca_fraction x {f}
end


create script_out /set_Cafraction_nmda
setfield /set_Cafraction_nmda command set_Cafraction



// make the control panel
make_control


// make the graph to display soma Vm and pass messages to the graph
make_Vmgraph

addmsg /cell/soma /data/voltage PLOT Vm *soma *black
addmsg /cell/trunk[6] /data/voltage PLOT Vm *trunk *red 


make_Gkgraph


/* comment out the line below to disable the cell display (faster)  */
make_xcell // create and display the xcell Vm

xcolorscale hot


//Granule cell stimulation
ce /

/* create presynaptic element for NMDA and AMPA channels */

create neutral presyn
setfield presyn z 0


//Stimulation patterns

//Spines on the tip of the terminal dendrites
addmsg /presyn /cell/periph11[20]/##[][TYPE=synchan] ACTIVATION z
addmsg /presyn /cell/periph12[17]/##[][TYPE=synchan] ACTIVATION z
addmsg /presyn /cell/periph21[5]/##[][TYPE=synchan] ACTIVATION z
addmsg /presyn /cell/periph22[8]/##[][TYPE=synchan] ACTIVATION z



//Setting maximum conductances
setfield /cell/##[][TYPE=compartment]/head[]/AMPA_NMDA gmax 1e-9
setfield /cell/##[][TYPE=compartment]/head[]/AMPA_NMDA/NMDA gmax 0.593e-9


//Ploting variables
addmsg /cell/periph12[17]/head[]/AMPA_NMDA/NMDA /data2/conductance PLOT Gk *GNMDA *green
addmsg /cell/periph12[17]/head[]/AMPA_NMDA /data2/conductance PLOT Gk *GAMPA *red


//Magnesium concentration
setfield /cell/##[][TYPE=Mg_block] CMg 1.2


reset









































 







