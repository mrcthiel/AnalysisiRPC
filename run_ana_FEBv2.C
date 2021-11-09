void run_ana_FEBv2(Int_t hv, Int_t sn, Int_t mt,/* Int_t muW1, Int_t muW2 ,*/ const char * loc){
            ana_FEBv2 t(hv, sn, mt,/* muW1, muW2,*/loc);
            t.Loop();
}	
