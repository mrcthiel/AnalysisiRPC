long int m_traw(long int m) {return m&0XFFFFFFFF;}
int m_fpga(long int m) {return (m>>32)&0X3;}
int m_channel(long int m) {return (m>>34)&0X3F;}
int m_side(long int m) {return (m>>40)&0X1;}
int m_strip(long int m) {return (m>>41)&0X3F;}
int m_pr_channel(long int m) {return (m>>47)&0X1F;}
float m_time(long int traw) {return traw*2.5/256;}
float m_tcor(long int traw, float tbc0) {return (traw-tbc0)*2.5/256;}

// Side 0 HR 1 LR
//int c_side(int ch) {return (ch>15)?ch%2:(ch+1)%2;}
int c_side(int ch) {return (ch)%2;}
int c_local_strip(int ch) {return (ch>15)?((ch-16)/2):16-(ch/2+1);}
int c_strip(int fpga, int ch) {return fpga*16+c_local_strip(ch);}
int c_petiroc(int ch) {return (ch>15)?(ch-16)*2:ch*2;}


// new FEBv2r2 maping
std::vector<int> m_strip_FEBv2r2(long int m) {

	int pr_ch  = m_pr_channel(m);
	int tdc_ch = m_channel(m);
	int fpga = m_fpga(m);

	std::vector<int> strip_side;
	strip_side.clear();

	if(fpga==0 && pr_ch==1 && tdc_ch==0){strip_side.push_back(15); strip_side.push_back(0);} 
	else if(fpga==0 && pr_ch==30 && tdc_ch==31){strip_side.push_back(15); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==2  && tdc_ch==1 ){strip_side.push_back(14); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==28 && tdc_ch==30){strip_side.push_back(14); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==4  && tdc_ch==2 ){strip_side.push_back(13); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==26 && tdc_ch==29){strip_side.push_back(13); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==6  && tdc_ch==3 ){strip_side.push_back(12); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==24 && tdc_ch==28){strip_side.push_back(12); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==8  && tdc_ch==4 ){strip_side.push_back(11); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==22 && tdc_ch==27){strip_side.push_back(11); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==10 && tdc_ch==5 ){strip_side.push_back(10); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==20 && tdc_ch==26){strip_side.push_back(10); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==12 && tdc_ch==6 ){strip_side.push_back(9); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==18 && tdc_ch==25){strip_side.push_back(9); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==14 && tdc_ch==7 ){strip_side.push_back(8); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==16 && tdc_ch==24){strip_side.push_back(8); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==16 && tdc_ch==8 ){strip_side.push_back(7); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==14 && tdc_ch==23){strip_side.push_back(7); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==18 && tdc_ch==9 ){strip_side.push_back(6); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==12 && tdc_ch==22){strip_side.push_back(6); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==20 && tdc_ch==10){strip_side.push_back(5); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==10 && tdc_ch==21){strip_side.push_back(5); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==22 && tdc_ch==11){strip_side.push_back(4); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==8  && tdc_ch==20){strip_side.push_back(4); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==24 && tdc_ch==12){strip_side.push_back(3); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==6  && tdc_ch==19){strip_side.push_back(3); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==26 && tdc_ch==13){strip_side.push_back(2); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==4  && tdc_ch==18){strip_side.push_back(2); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==28 && tdc_ch==14){strip_side.push_back(1); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==2  && tdc_ch==17){strip_side.push_back(1); strip_side.push_back(1);}
	else if(fpga==0 && pr_ch==30 && tdc_ch==15){strip_side.push_back(0); strip_side.push_back(0);}
	else if(fpga==0 && pr_ch==1  && tdc_ch==16){strip_side.push_back(0); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==1  && tdc_ch==0 ){strip_side.push_back(31); strip_side.push_back(0);} 
	else if(fpga==1 && pr_ch==30 && tdc_ch==31){strip_side.push_back(31); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==2  && tdc_ch==1 ){strip_side.push_back(30); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==28 && tdc_ch==30){strip_side.push_back(30); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==4  && tdc_ch==2 ){strip_side.push_back(29); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==26 && tdc_ch==29){strip_side.push_back(29); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==6  && tdc_ch==3 ){strip_side.push_back(28); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==24 && tdc_ch==28){strip_side.push_back(28); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==8  && tdc_ch==4 ){strip_side.push_back(27); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==22 && tdc_ch==27){strip_side.push_back(27); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==10 && tdc_ch==5 ){strip_side.push_back(26); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==20 && tdc_ch==26){strip_side.push_back(26); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==12 && tdc_ch==6 ){strip_side.push_back(25); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==18 && tdc_ch==25){strip_side.push_back(25); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==14 && tdc_ch==7 ){strip_side.push_back(24); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==16 && tdc_ch==24){strip_side.push_back(24); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==16 && tdc_ch==8 ){strip_side.push_back(23); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==14 && tdc_ch==23){strip_side.push_back(23); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==18 && tdc_ch==9 ){strip_side.push_back(22); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==12 && tdc_ch==22){strip_side.push_back(22); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==20 && tdc_ch==10){strip_side.push_back(21); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==10 && tdc_ch==21){strip_side.push_back(21); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==22 && tdc_ch==11){strip_side.push_back(20); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==8  && tdc_ch==20){strip_side.push_back(20); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==24 && tdc_ch==12){strip_side.push_back(19); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==6  && tdc_ch==19){strip_side.push_back(19); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==26 && tdc_ch==13){strip_side.push_back(18); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==4  && tdc_ch==18){strip_side.push_back(18); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==28 && tdc_ch==14){strip_side.push_back(17); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==2  && tdc_ch==17){strip_side.push_back(17); strip_side.push_back(1);}
	else if(fpga==1 && pr_ch==30 && tdc_ch==15){strip_side.push_back(16); strip_side.push_back(0);}
	else if(fpga==1 && pr_ch==1  && tdc_ch==16){strip_side.push_back(16); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==1  && tdc_ch==0 ){strip_side.push_back(47); strip_side.push_back(0);} 
	else if(fpga==2 && pr_ch==30 && tdc_ch==31){strip_side.push_back(47); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==2  && tdc_ch==1 ){strip_side.push_back(46); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==28 && tdc_ch==30){strip_side.push_back(46); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==4  && tdc_ch==2 ){strip_side.push_back(45); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==26 && tdc_ch==29){strip_side.push_back(45); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==6  && tdc_ch==3 ){strip_side.push_back(44); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==24 && tdc_ch==28){strip_side.push_back(44); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==8  && tdc_ch==4 ){strip_side.push_back(43); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==22 && tdc_ch==27){strip_side.push_back(43); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==10 && tdc_ch==5 ){strip_side.push_back(42); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==20 && tdc_ch==26){strip_side.push_back(42); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==12 && tdc_ch==6 ){strip_side.push_back(41); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==18 && tdc_ch==25){strip_side.push_back(41); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==14 && tdc_ch==7 ){strip_side.push_back(40); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==16 && tdc_ch==24){strip_side.push_back(40); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==16 && tdc_ch==8 ){strip_side.push_back(39); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==14 && tdc_ch==23){strip_side.push_back(39); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==18 && tdc_ch==9 ){strip_side.push_back(38); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==12 && tdc_ch==22){strip_side.push_back(38); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==20 && tdc_ch==10){strip_side.push_back(37); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==10 && tdc_ch==21){strip_side.push_back(37); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==22 && tdc_ch==11){strip_side.push_back(36); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==8  && tdc_ch==20){strip_side.push_back(36); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==24 && tdc_ch==12){strip_side.push_back(35); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==6  && tdc_ch==19){strip_side.push_back(35); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==26 && tdc_ch==13){strip_side.push_back(34); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==4  && tdc_ch==18){strip_side.push_back(34); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==28 && tdc_ch==14){strip_side.push_back(33); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==2  && tdc_ch==17){strip_side.push_back(33); strip_side.push_back(1);}
	else if(fpga==2 && pr_ch==30 && tdc_ch==15){strip_side.push_back(32); strip_side.push_back(0);}
	else if(fpga==2 && pr_ch==1  && tdc_ch==16){strip_side.push_back(32); strip_side.push_back(1);}

	return strip_side;
}
