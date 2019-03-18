//genesis

// CONSTANTS


//genesis
/* Modified from Traub91_proto*/

/* Constants*/
    float PI = 3.14159
    float CM = 0.01
    float RM = 1.20
    float RA = 2.00
    float EREST_ACT = -0.070

float ECA = 0.140 + EREST_ACT 
float SOMA_A = 3.320e-9       // soma area in square meters

    /* Sizes */
    /* Spine*/
    float head_dia = 0.8e-6 //on average to Granule cells peduculated spines (Woolf, Shepherd, Greer, 1991)
    float head_len = 0.8e-6 //on average to Granule cells peduculated spines (Woolf, Shepherd, Greer, 1991)
    float neck_dia = 0.23e-6 //on average to Granule cells peduculated spines (Woolf, Shepherd, Greer, 1991)
    float neck_len = 1.9e-6 //on average to Granule cells peduculated spines (Woolf, Shepherd, Greer, 1991)

    /* soma */
    float soma_dia = 10.0e-6
    /* Shell thickness */
    float thick = 0.1e-6


    /* Synaptic channels */
//Parameters from Davison, Feng e Brown, J.Neurophysiol., vol90, pg.1921-1935, 2003.
//From spines of granule cells of olfactory bulb
    float G_NMDA = 0.593e-9	// maximum conductance
    float E_NMDA = 0.0		// reversal potential
    float tau1_NMDA = 52.0e-3	// open time constant
    float tau2_NMDA = 343.0e-3	// close time constant
    float CMg = 1.2	// Magnesium concentration for magnesium block
    float eta = 0.2801
    float gamma = 62
    float fraction = 0.185	// relative fraction of Ca current flowing into shell

    float G_AMPA = 1.0e-9	// maximum conductance
    float E_AMPA = 0.0		// reversal potential
    float tau1_AMPA = 2.0e-3	// open time constant
    float tau2_AMPA = 5.5e-3	// close time constant
 
    float G_GABA_A= 1.1e-9	// maximum conductance
    float E_GABA_A= -0.075 // reversal potential
    float tau1_GABA_A= 4e-3// open time constant
    float tau2_GABA_A= 18e-3// close time constant



//========================================================================
//                      AMPA Channel, NMDA channel with Mg_block
//========================================================================


function make_AMPA_NMDA
    /* add a non-NMDA channel: is always activated together with the NMDA
    ** channel */
    create synchan AMPA_NMDA  //cria o AMPA
    setfield AMPA_NMDA gmax {G_AMPA} Ek {E_AMPA} tau1 {tau1_AMPA} tau2 {tau2_AMPA}

// NMDA Channeland Mg_block

ce AMPA_NMDA
    /* add a NMDA channel: is used to compute channel conductance only */
    create synchan NMDA  //cria o NMDA
    setfield NMDA gmax {G_NMDA} Ek {E_NMDA} tau1 {tau1_NMDA} tau2 {tau2_NMDA}

    /* add the Mg block: the blocked NMDA current is used to compute voltage */
    create Mg_block Mg_block
    setfield Mg_block CMg {CMg}  Ek {E_NMDA} Zk 2 \
            KMg_A {1/eta} \ \\ *({exp {EREST_ACT*gamma}})} \
            KMg_B {1.0/gamma}

   addmsg NMDA Mg_block CHANNEL Gk Ek
	
	//Ca fraction
	create neutral Ca_fraction
	setfield Ca_fraction x 0.185
	
	//iCa fraction
	create calculator ICa_fraction

ce ..

   addfield AMPA_NMDA addmsg1
   setfield AMPA_NMDA addmsg1 ".. ./Mg_block VOLTAGE Vm"

   addfield AMPA_NMDA addmsg2
   setfield AMPA_NMDA addmsg2 "./Mg_block .. CHANNEL Gk Ek"

   addfield AMPA_NMDA addmsg3
   setfield AMPA_NMDA addmsg3 "./Mg_block ./ICa_fraction SUM Ik" 

   addfield AMPA_NMDA addmsg4
   setfield AMPA_NMDA addmsg4 "./Ca_fraction ./ICa_fraction MULTIPLY x" 
   
end


//========================================================================
//                      GABA_A Channel
//========================================================================


function make_GABA_A
    create synchan GABA_A  //cria o GABA_A
    setfield GABA_A gmax {G_GABA_A} Ek {E_GABA_A} tau1 {tau1_GABA_A} tau2 {tau2_GABA_A}
end



//========================================================================
//                      Ca conc 
//========================================================================

function make_Ca_conc
        if ({exists Ca_conc})
                return
        end
        
        create Ca_concen Ca_conc
        setfield Ca_conc \
			tau     0.0011 \ // sec (Egger and Stroh, 2009)
                B       26e11 \  // Curr to conc for soma
                Ca_base 50e-6    //0.05 uM (Egger and Stroh, 2009)


        addfield Ca_conc addmsg1
        setfield Ca_conc addmsg1 "../AMPA_NMDA/ICa_fraction . I_Ca output"


end






//========================================================================
//                      Tabulated Ca ChannelTraub 91
//========================================================================

function make_Ca
        if ({exists Ca})
                return
        end

        create  tabchannel      Ca
                setfield        ^       \
                Ek              {ECA}   \               //      V
                Gbar            { 40 * SOMA_A }      \  //      S
                Ik              0       \               //      A
                Gk              0       \               //      S
                Xpower  2       \
                Ypower  1       \
                Zpower  0


        setupalpha Ca X 1.6e3  \
                 0 1.0 {-1.0 * (0.065 + EREST_ACT) } -0.01389       \
                 {-20e3 * (0.0511 + EREST_ACT) }  \
                 20e3 -1.0 {-1.0 * (0.0511 + EREST_ACT) } 5.0e-3 


        float   xmin = -0.1
        float   xmax = 0.05
        int     xdivs = 49
	call Ca TABCREATE Y {xdivs} {xmin} {xmax}

        int i
        float x,dx,y
        dx = (xmax - xmin)/xdivs
        x = xmin
        for (i = 0 ; i <= {xdivs} ; i = i + 1)
	    if (x > EREST_ACT)
                y = 5.0*{exp {-50*(x - EREST_ACT)}}
	    else
		y = 5.0
	    end
            setfield Ca Y_A->table[{i}] {y}
            setfield Ca Y_B->table[{i}] 5.0
            x = x + dx
        end

           setfield Ca Y_A->calc_mode 0   Y_B->calc_mode 0
           call Ca TABFILL Y 3000 0
end



















