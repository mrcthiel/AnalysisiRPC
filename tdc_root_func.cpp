long int m_traw(long int m) {return m&0XFFFFFFFF;}
int m_fpga(long int m) {return (m>>32)&0X3;}
int m_channel(long int m) {return (m>>34)&0X3F;}
int m_side(long int m) {return (m>>40)&0X1;}
int m_strip(long int m) {return (m>>41)&0X3F;}
int m_pr_channel(long int m) {return (m>>47)&0X1F;}
float m_time(long int traw) {return traw*2.5/256;}
float m_tcor(long int traw, float tbc0) {return (traw-tbc0)*2.5/256;}

// Side 0 LR 1 HR
int c_side(int ch) {return (ch>15)?ch%2:(ch+1)%2;}
int c_local_strip(int ch) {return (ch>15)?((ch-16)/2):16-(ch/2+1);}
int c_strip(int fpga, int ch) {return fpga*16+c_local_strip(ch);}
int c_petiroc(int ch) {return (ch>15)?(ch-16)*2:ch*2;}

