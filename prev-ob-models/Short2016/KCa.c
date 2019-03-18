/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 
#define _threadargsprotocomma_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkbar _p[0]
#define ik _p[1]
#define Yvdep _p[2]
#define Yconcdep _p[3]
#define Y _p[4]
#define ek _p[5]
#define cai _p[6]
#define DY _p[7]
#define _g _p[8]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_concdep(void);
 static void _hoc_rate(void);
 static void _hoc_vdep(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_Ikca", _hoc_setdata,
 "concdep_Ikca", _hoc_concdep,
 "rate_Ikca", _hoc_rate,
 "vdep_Ikca", _hoc_vdep,
 0, 0
};
 /* declare global and static user variables */
#define Ybeta Ybeta_Ikca
 double Ybeta = 0.05;
#define Yalpha Yalpha_Ikca
 double Yalpha = 0;
#define usetable usetable_Ikca
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gkbar_Ikca", 0, 1e+009,
 "usetable_Ikca", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Ybeta_Ikca", "/ms",
 "Yalpha_Ikca", "/ms",
 "gkbar_Ikca", "mho/cm2",
 "ik_Ikca", "mA/cm2",
 "Yconcdep_Ikca", "/ms",
 0,0
};
 static double Y0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Ybeta_Ikca", &Ybeta_Ikca,
 "Yalpha_Ikca", &Yalpha_Ikca,
 "usetable_Ikca", &usetable_Ikca,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"Ikca",
 "gkbar_Ikca",
 0,
 "ik_Ikca",
 "Yvdep_Ikca",
 "Yconcdep_Ikca",
 0,
 "Y_Ikca",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.12;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _KCa_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 9, 5);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Ikca e:/Dropbox/Research/Thesis/MODELS/Short2016/KCa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_Yvdep;
 static double *_t_Yconcdep;
static int _reset;
static char *modelname = "Calcium dependent potassium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_concdep(double);
static int _f_vdep(double);
static int concdep(double);
static int rate(double, double);
static int vdep(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_concdep(double);
 static void _n_vdep(double);
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rate ( _threadargscomma_ v , cai ) ;
   DY = Yalpha * ( 1.0 - Y ) - Ybeta * Y ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rate ( _threadargscomma_ v , cai ) ;
 DY = DY  / (1. - dt*( (Yalpha)*(( ( - 1.0 ) )) - (Ybeta)*(1.0) )) ;
 return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
   rate ( _threadargscomma_ v , cai ) ;
    Y = Y + (1. - exp(dt*((Yalpha)*(( ( - 1.0 ) )) - (Ybeta)*(1.0))))*(- ( (Yalpha)*(( 1.0 )) ) / ( (Yalpha)*(( ( - 1.0) )) - (Ybeta)*(1.0) ) - Y) ;
   }
  return 0;
}
 
static int  rate (  double _lv , double _lcai ) {
   vdep ( _threadargscomma_ _lv ) ;
   concdep ( _threadargscomma_ _lcai ) ;
   Yalpha = Yvdep * Yconcdep ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   _r = 1.;
 rate (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 static double _mfac_vdep, _tmin_vdep;
 static void _check_vdep();
 static void _check_vdep() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_vdep =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_vdep)/100.; _mfac_vdep = 1./_dx;
   for (_i=0, _x=_tmin_vdep; _i < 101; _x += _dx, _i++) {
    _f_vdep(_x);
    _t_Yvdep[_i] = Yvdep;
   }
  }
 }

 static int vdep(double _lv){ _check_vdep();
 _n_vdep(_lv);
 return 0;
 }

 static void _n_vdep(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_vdep(_lv); return; 
}
 _xi = _mfac_vdep * (_lv - _tmin_vdep);
 if (isnan(_xi)) {
  Yvdep = _xi;
  return;
 }
 if (_xi <= 0.) {
 Yvdep = _t_Yvdep[0];
 return; }
 if (_xi >= 100.) {
 Yvdep = _t_Yvdep[100];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 Yvdep = _t_Yvdep[_i] + _theta*(_t_Yvdep[_i+1] - _t_Yvdep[_i]);
 }

 
static int  _f_vdep (  double _lv ) {
   Yvdep = exp ( ( _lv * 1.0 - 65.0 ) / 27.0 ) ;
    return 0; }
 
static void _hoc_vdep(void) {
  double _r;
    _r = 1.;
 vdep (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_concdep, _tmin_concdep;
 static void _check_concdep();
 static void _check_concdep() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_concdep =  0.0 ;
   _tmax =  50.0 ;
   _dx = (_tmax - _tmin_concdep)/2000.; _mfac_concdep = 1./_dx;
   for (_i=0, _x=_tmin_concdep; _i < 2001; _x += _dx, _i++) {
    _f_concdep(_x);
    _t_Yconcdep[_i] = Yconcdep;
   }
  }
 }

 static int concdep(double _lcai){ _check_concdep();
 _n_concdep(_lcai);
 return 0;
 }

 static void _n_concdep(double _lcai){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_concdep(_lcai); return; 
}
 _xi = _mfac_concdep * (_lcai - _tmin_concdep);
 if (isnan(_xi)) {
  Yconcdep = _xi;
  return;
 }
 if (_xi <= 0.) {
 Yconcdep = _t_Yconcdep[0];
 return; }
 if (_xi >= 2000.) {
 Yconcdep = _t_Yconcdep[2000];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 Yconcdep = _t_Yconcdep[_i] + _theta*(_t_Yconcdep[_i+1] - _t_Yconcdep[_i]);
 }

 
static int  _f_concdep (  double _lcai ) {
   if ( fabs ( _lcai - 0.015 ) < 1e-5 ) {
     Yconcdep = 500.0 / ( 1.0 / 0.0013 + 0.5 * ( 0.015 - _lcai ) / pow( 0.0013 , 2.0 ) ) ;
     }
   else {
     Yconcdep = 500.0 * ( 0.015 - _lcai * 1.0 ) / ( exp ( ( 0.015 - _lcai * 1.0 ) / 0.0013 ) - 1.0 ) ;
     }
    return 0; }
 
static void _hoc_concdep(void) {
  double _r;
    _r = 1.;
 concdep (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  Y = Y0;
 {
   rate ( _threadargscomma_ v , cai ) ;
   Y = Yalpha / ( Yalpha + Ybeta ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gkbar * Y * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ek = _ion_ek;
  cai = _ion_cai;
 { {
 for (; t < _break; t += dt) {
 error =  state();
 if(error){fprintf(stderr,"at line 55 in file KCa.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(Y) - _p;  _dlist1[0] = &(DY) - _p;
   _t_Yvdep = makevector(101*sizeof(double));
   _t_Yconcdep = makevector(2001*sizeof(double));
_first = 0;
}
